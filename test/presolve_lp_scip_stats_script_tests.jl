using Test
import MathOptInterface as MOI
import QIPresolve
import SCIP

include(joinpath(@__DIR__, "..", "scripts", "presolve_lp_scip_stats.jl"))
const LPStatsScript = Main.PresolveLpScipStatsScript

const LP_STATS_FEASIBLE = """
Maximize
 obj: x1 + x2
Subject To
 c1: x1 + x2 <= 1
Binary
 x1
 x2
End
"""

function check_lp_stats_timings(result)
    times = (result.wall_time_sec, result.qip_presolve_time_sec,
        result.scip_presolve_time_sec, result.scip_solving_time_sec)
    @test all(t -> isfinite(t) && t >= 0.0, times)
    @test result.wall_time_sec >= result.qip_presolve_time_sec
    # Allow timer resolution differences between Julia and SCIP.
    @test result.wall_time_sec + 0.01 >= sum(times[2:4])
end

@testset "LP SCIP comparison prints only the requested metrics" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "tiny.lp")
        write(lp_path, LP_STATS_FEASIBLE)
        output = IOBuffer()
        results = LPStatsScript.main([lp_path]; io = output)
        blocks = split(strip(String(take!(output))), "\n\n")

        @test length(blocks) == 2
        labels = ("SCIP only", "Parity/residue + SCIP")
        metric_labels = ["Full wall time (s)", "Parity/residue presolve time (s)",
            "SCIP presolve time (s)", "SCIP solving time (s)", "SCIP status"]
        for (block, label, result) in zip(blocks, labels, results)
            lines = split(block, '\n')
            @test length(lines) == 6
            @test lines[1] == label
            @test [strip(first(split(line, ':'; limit = 2))) for line in lines[2:end]] == metric_labels
            @test result.status == "SCIP_STATUS_OPTIMAL"
            @test last(lines) == "  SCIP status: SCIP_STATUS_OPTIMAL"
            check_lp_stats_timings(result)
        end
        @test results.original.qip_presolve_time_sec == 0.0
        @test readdir(dir) == ["tiny.lp"]
        @test read(lp_path, String) == LP_STATS_FEASIBLE
    end
end

@testset "LP SCIP comparison skips SCIP when QIPresolve proves infeasibility" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "infeasible.lp")
        write(lp_path, replace(LP_STATS_FEASIBLE, "x1 + x2 <= 1" => "2 x1 + 2 x2 = 1"))
        output = IOBuffer()
        results = LPStatsScript.main([lp_path]; io = output)
        @test results.original.status == "SCIP_STATUS_INFEASIBLE"
        @test results.presolved.status == "INFEASIBLE (QIPresolve; SCIP skipped)"
        @test results.presolved.scip_presolve_time_sec == 0.0
        @test results.presolved.scip_solving_time_sec == 0.0
        @test occursin("SCIP status: INFEASIBLE (QIPresolve; SCIP skipped)", String(take!(output)))
        check_lp_stats_timings(results.original)
        check_lp_stats_timings(results.presolved)
    end
end

@testset "LP SCIP comparison solves models completely reduced by presolve" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "reduced.lp")
        write(lp_path, replace(LP_STATS_FEASIBLE, "x1 + x2 <= 1" => "x1 + 2 x2 = 1"))
        lp_model = LPStatsScript.read_lp(lp_path)
        model = LPStatsScript.build_core_model(lp_model)
        @test all(v -> v.lb == 0.0 && v.ub == 1.0, values(model.vars))
        QIPresolve.presolve!(model)
        @test isempty(model.vars)
        @test isempty(model.cons)
        @test !model.infeasible

        results = LPStatsScript.main([lp_path]; io = IOBuffer())
        @test results.original.status == results.presolved.status == "SCIP_STATUS_OPTIMAL"
        check_lp_stats_timings(results.presolved)
    end
end

@testset "LP SCIP comparison supports quadratic constraints" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "quadratic.lp")
        write(lp_path, replace(LP_STATS_FEASIBLE, "x1 + x2 <= 1" => "[ 2 x1 * x2 ] <= 1"))
        results = LPStatsScript.main([lp_path]; io = IOBuffer())
        @test results.original.status == results.presolved.status == "SCIP_STATUS_OPTIMAL"
        check_lp_stats_timings(results.presolved)
    end
end

@testset "SCIP time limit and timing counters" begin
    optimizer = SCIP.Optimizer()
    try
        LPStatsScript.configure_scip!(optimizer)
        @test MOI.get(optimizer, MOI.TimeLimitSec()) == 1800.0
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("timing/clocktype")) == 2
        @test MOI.get(optimizer, MOI.Silent())
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("presolving/maxrounds")) == -1
        x = MOI.add_variable(optimizer)
        MOI.add_constraint(optimizer, x, MOI.ZeroOne())
        @test LPStatsScript.scip_metrics(optimizer).scip_presolve_time_sec == 0.0

        # Override the limit only in this test to exercise an immediate stop.
        MOI.set(optimizer, MOI.TimeLimitSec(), 0.0)
        MOI.optimize!(optimizer)
        metrics = LPStatsScript.scip_metrics(optimizer)
        @test metrics.status == "SCIP_STATUS_TIMELIMIT"
        # SCIP may briefly enter presolving before checking even a zero limit.
        @test isfinite(metrics.scip_presolve_time_sec) && metrics.scip_presolve_time_sec >= 0.0
        @test isfinite(metrics.scip_solving_time_sec)
        @test metrics.scip_solving_time_sec >= 0.0
    finally
        SCIP.free_scip(optimizer.inner)
    end

    optimizer = SCIP.Optimizer()
    try
        LPStatsScript.configure_scip!(optimizer)
        x = MOI.add_variable(optimizer)
        MOI.add_constraint(optimizer, x, MOI.ZeroOne())
        MOI.optimize!(optimizer)
        metrics = LPStatsScript.scip_metrics(optimizer)
        @test metrics.status == "SCIP_STATUS_OPTIMAL"
        @test metrics.scip_presolve_time_sec == SCIP.SCIPgetPresolvingTime(optimizer.inner)
        @test metrics.scip_presolve_time_sec + metrics.scip_solving_time_sec ≈
            SCIP.SCIPgetSolvingTime(optimizer.inner)
    finally
        SCIP.free_scip(optimizer.inner)
    end
end

@testset "LP SCIP comparison CLI and unsupported inputs" begin
    for flag in ("-h", "--help")
        output = IOBuffer()
        @test LPStatsScript.main([flag]; io = output) === nothing
        @test occursin("Usage:", String(take!(output)))
    end
    for args in (String[], ["one.lp", "two.lp"], ["--silent"],
            ["one.lp", "--output", "stats.json"], ["--help", "one.lp"])
        @test_throws ErrorException LPStatsScript.main(args; io = IOBuffer())
    end
    mktempdir() do dir
        @test_throws "LP file not found" LPStatsScript.main([joinpath(dir, "missing.lp")])
        lp_path = joinpath(dir, "unsupported.lp")
        write(lp_path, replace(LP_STATS_FEASIBLE, "Binary\n x1\n x2\n" => ""))
        @test_throws "only integer and binary variables" LPStatsScript.main([lp_path])
        write(lp_path, replace(LP_STATS_FEASIBLE, "obj: x1 + x2" => "obj: [ 2 x1 ^ 2 ] / 2"))
        @test_throws "affine objectives only" LPStatsScript.main([lp_path])
    end
end
