module MainBenchmarkComparisonTests
using Test
using CSV
import SCIP
import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC
include(joinpath(@__DIR__, "..", "scripts", "presolve_main_benchmark_comparison.jl"))
const C = MainBenchmarkPresolveComparison
const QT = Tuple{Float64, Int, Int}
const LT = Tuple{Float64, Int}
empty_expr() = PC.QuadExpr(QT[], LT[])
model(cons; ub = 10.0) = PC.QPModel(Dict(i => PC.IntVar(0.0, ub) for i in 1:2), cons, empty_expr(), :min)
row(id, lhs, rhs; coeffs = [(2.0, 1), (3.0, 2)]) = PC.Constraint(id, PC.QuadExpr(QT[], coeffs), lhs, rhs)

@testset "comparison metric definitions and provenance" begin
    @test C.domain_reduction(log(100), log(10)) ≈ 0.5
    @test C.domain_reduction(0.0, 0.0) == 0.0
    @test C.domain_reduction(log(100), 0.0) == 1.0
    @test C.log_domain_sum(model(PC.Constraint[]; ub = 9.0)) ≈ log(100)
    @test C.average_contributions(Float64[]) == 0.0
    for (lhs, rhs, width, expected) in ((1.0, 5.0, 2.0, 0.5),
            (1.0, 5.0, -1.0, 1.25), (1.0, 1.0, 0.0, 0.0), (-Inf, 5.0, 0.0, 0.0))
        original = row(1, lhs, rhs)
        baseline = only(C.bound_baselines(model([original])))
        @test C.bound_contribution(baseline, width) == expected
        @test C.bound_contribution(baseline, width; removed = true) == 1.0
    end
    original = row(1, 0.0, 12.0; coeffs = [(2.0, 1), (4.0, 2)])
    m = model([original])
    bounds = C.bound_baselines(m)
    PC.scale_gcd!(original)
    @test C.core_contributions(m, bounds) == [0.0]
    original.rhs -= 1.0
    @test only(C.core_contributions(m, bounds)) ≈ 1 / 6
    PC.affine_transform!(original, 1, 1.0, 2.0)
    @test only(C.core_contributions(m, bounds)) ≈ 1 / 6
    PC._scale_constraint_by_two!(original)
    @test only(C.core_contributions(m, bounds)) ≈ 1 / 6
    # Identity, rather than id, distinguishes removed rows from new helper rows.
    m.cons = [row(1, -10.0, 10.0)]
    @test C.core_contributions(m, bounds) == [1.0]

    m = model([row(1, 0.0, 10.0), row(2, 2.0, 8.0)])
    bounds = C.bound_baselines(m)
    PC.aggregate_parallel_constraints!(m)
    @test C.core_contributions(m, bounds) ≈ [0.4, 1.0]
    split = model([row(1, -Inf, 8.0), row(2, 2.0, Inf), row(3, 3.0, 3.0)])
    C.rejoin_ranges!(split)
    @test length(split.cons) == 2
    @test any(c -> c.lhs == 2.0 && c.rhs == 8.0, split.cons)
    split_neg = model([row(1, -Inf, 8.0), row(2, -Inf, -2.0; coeffs = [(-2.0, 1), (-3.0, 2)])])
    C.rejoin_ranges!(split_neg)
    @test length(split_neg.cons) == 1
    @test only(split_neg.cons).rhs - only(split_neg.cons).lhs == 6.0
end

@testset "SCIP comparison and strategy independence" begin
    config = C.Config()
    m = model([row(1, 1.0, 4.0)]; ub = 1.0)
    before = PC._model_state_signature(m)
    iso = C.run_core(m, config; enable_parity = false, enable_residue = true)
    @test only(iso.contributions) ≈ 2 / 3
    @test iso.status == "reduced"
    @test PC._model_state_signature(m) == before
    scip = C.run_scip(m, C.bound_baselines(m), config; diagnostics = devnull)
    @test scip.nodes == 0
    @test scip.status == "feasible"
    # SCIP can solve a feasibility objective before deleting the row. A solved
    # status alone must not turn an unchanged interval into complete tightening.
    @test scip.contributions == [0.0]
    @test PC._model_state_signature(m) == before

    infeasible = model([row(1, 1.0, 1.0), row(2, 2.0, 2.0)]; ub = 1.0)
    scip = C.run_scip(infeasible, C.bound_baselines(infeasible), config; diagnostics = devnull)
    @test scip.status == "infeasible"
    @test scip.log_domain == 0.0
    @test scip.nodes == 0

    # Full-stage accounting keeps removals from QIPresolve in its original denominator.
    partial = model([row(1, 1.0, 4.0), row(2, 1.0, 1.0)]; ub = 1.0)
    baselines = C.bound_baselines(partial)
    pop!(partial.cons)
    full = C.run_scip(partial, baselines, config; diagnostics = devnull)
    @test full.contributions == [0.0, 1.0]
    eliminated = model([row(1, 3.0, 11.0; coeffs = [(2.0, 1), (4.0, 2)])])
    removed = C.run_scip(eliminated, C.bound_baselines(eliminated), config; diagnostics = devnull)
    @test removed.status == "feasible"
    @test removed.contributions == [1.0]

    p = C.Polynomial((UInt(1), UInt(0)) => 2.0, C.CONSTANT => 3.0)
    q = C.Polynomial((UInt(1), UInt(0)) => -4.0, C.CONSTANT => 10.0)
    @test C.proportional_scale(p, q) == -2.0
    q[(UInt(2), UInt(0))] = 1.0
    @test C.proportional_scale(p, q) === nothing
end

@testset "SCIP replacements retain constraint provenance" begin
    mktempdir() do dir
        settings = joinpath(dir, "trace.set")
        optimizer = SCIP.Optimizer()
        try
            SCIP.@SCIP_CALL SCIP.SCIPsetHeuristics(optimizer.inner, SCIP.SCIP_PARAMSETTING_OFF, 1)
            SCIP.@SCIP_CALL SCIP.SCIPsetBoolParam(optimizer.inner, "misc/allowstrongdualreds", 0)
            SCIP.@SCIP_CALL SCIP.SCIPsetBoolParam(optimizer.inner, "misc/allowweakdualreds", 0)
            SCIP.@SCIP_CALL SCIP.SCIPwriteParams(optimizer.inner, settings, 0, 1)
        finally
            SCIP.free_scip(optimizer.inner)
        end
        config = C.Config(scip_config = settings)
        cases = [(model([row(1, 3.0, 11.0; coeffs = [(2.0, 1), (4.0, 2)])]), "varbound", 0.25),
            (PC.QPModel(Dict(i => PC.IntVar(0.0, 1.0) for i in 1:3),
                [row(1, 0.0, 1.0; coeffs = [(1.0, i) for i in 1:3])], empty_expr(), :min), "setppc", 0.0)]
        for (m, expected_handler, tightening) in cases
            optimizer, _ = C.build_native_model(m)
            try
                SCIP.@SCIP_CALL SCIP.SCIPreadParams(optimizer.inner, settings)
                C.MOI.set(optimizer, C.MOI.Silent(), true)
                SCIP.@SCIP_CALL SCIP.SCIPpresolve(optimizer.inner)
                conss = unsafe_wrap(Array, SCIP.SCIPgetConss(optimizer.inner), SCIP.SCIPgetNConss(optimizer.inner))
                @test length(conss) == 1
                @test unsafe_string(SCIP.SCIPconshdlrGetName(SCIP.SCIPconsGetHdlr(only(conss)))) == expected_handler
            finally
                SCIP.free_scip(optimizer.inner)
            end
            result = C.run_scip(m, C.bound_baselines(m), config; diagnostics = devnull)
            @test result.status == "reduced"
            @test result.nodes == 0
            @test only(result.contributions) ≈ tightening
        end
        # Treat the untracked duplicate as an input helper. It must not be
        # mistaken for a newly created successor when SCIP removes the original.
        duplicates = model([row(i, 3.0, 11.0; coeffs = [(2.0, 1), (4.0, 2)]) for i in 1:2])
        contributions = [only(C.run_scip(duplicates, [baseline], config;
            diagnostics = devnull).contributions) for baseline in C.bound_baselines(duplicates)]
        @test sort(contributions) ≈ [0.25, 1.0]
    end
end

@testset "benchmark CSV integration" begin
    mktempdir() do dir
        input = joinpath(dir, "input")
        output = joinpath(dir, "output")
        mkpath(input)
        write(joinpath(input, "random_instances.csv"),
            "instance_name,file_name,subtype\nrandom_unit_001,random_unit_001.lp,unit\nrandom_unit_002,random_unit_002.lp,unit\n")
        for (index, rhs) in ((1, 1), (2, 4))
            write(joinpath(input, "random_unit_00$index.lp"), """
            minimize
            obj:
            subject to
            c: x + 2 y = $rhs
            bounds
            0 <= x <= 1
            0 <= y <= 1
            general
            x y
            end
            """)
        end
        config = C.parse_args(["--input-dir", input, "--output-dir=$output"])
        results = C.run_experiment(config)
        @test length(results) == 2
        @test results[1].status_comb == results[1].status_full == "feasible"
        @test results[2].status_comb == results[2].status_full == "infeasible"
        @test results[1].bound_tightening_comb == results[1].bound_tightening_full == 1.0
        per_instance = collect(CSV.File(joinpath(output, "per_instance.csv")))
        aggregated = collect(CSV.File(joinpath(output, "aggregated.csv")))
        @test length(per_instance) == 2
        @test length(aggregated) == 1
        @test aggregated[1].instance_count == 2
        @test !hasproperty(aggregated[1], :status_comb)
        @test all(hasproperty(per_instance[1], key) for key in C.METRICS)
        @test aggregated[1].bound_tightening_full ≈ sum(r.bound_tightening_full for r in results) / 2
        @test C.parse_args(["--limit=1"]).limit == 1
        @test_throws ErrorException C.parse_args(["--limit=0"])
        @test_throws ErrorException C.parse_args(["--unknown=1"])
    end
end
end # module
