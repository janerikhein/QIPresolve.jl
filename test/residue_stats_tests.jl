using Test
import QIPresolve
import QIPresolve.PresolvingCore as PC

const ResidueStatsQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const ResidueStatsLinTerm = Tuple{Float64, PC.VarId}
const RESIDUE_STATS_NEXT_CON_ID = Ref(0)

residue_stats_next_con_id() = (RESIDUE_STATS_NEXT_CON_ID[] += 1)
residue_stats_empty_objective() = PC.QuadExpr(ResidueStatsQuadTerm[], ResidueStatsLinTerm[])

function residue_stats_model(vars, cons)
    return PC.QPModel(vars, cons, residue_stats_empty_objective(), :min)
end

function residue_stats_constraint(quad_terms, lin_terms, lhs::Float64, rhs::Float64)
    return PC.Constraint(
        residue_stats_next_con_id(),
        PC.QuadExpr(ResidueStatsQuadTerm[quad_terms...], ResidueStatsLinTerm[lin_terms...]),
        lhs,
        rhs,
    )
end

@testset "ResidueStats defaults and visibility" begin
    defaults = PC.ResidueStats()

    for field in fieldnames(PC.ResidueStats)
        value = getfield(defaults, field)
        if value isa Integer
            @test value == 0
        end
    end
    @test defaults.residue_presolve_time == 0.0
    @test defaults.avg_relative_bound_tightening == 0.0
    @test !defaults.infeasibility_detected

    @test isdefined(PC, :ResidueStats)
    @test !isdefined(QIPresolve, :ResidueStats)
end

@testset "ResidueStats records single constraint bound tightening" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(1.0, 4.0))
    con = residue_stats_constraint([], [(2.0, 1)], 1.0, 9.0)
    model = residue_stats_model(vars, [con])

    result = PC.residue_presolve_constraint!(model, con, [2])
    stats = result.residue_stats

    @test result.changed
    @test result.processed
    @test !result.tightened_to_equality
    @test !result.infeasible
    @test con.lhs == 2.0
    @test con.rhs == 8.0
    @test stats.num_constraints_processed == 1
    @test stats.num_constraints_tightened == 1
    @test stats.num_constraints_eq_tightened == 0
    @test stats.num_lower_bounds_tightened == 1
    @test stats.num_upper_bounds_tightened == 1
    @test stats.num_moduli_evaluated == 1
    @test stats.num_useful_moduli == 1
    @test stats.avg_relative_bound_tightening ≈ 0.25
    @test stats.residue_presolve_time >= 0.0
    @test !stats.infeasibility_detected
end

@testset "ResidueStats records equality tightening" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 2.0))
    con = residue_stats_constraint([(2.0, 1, 1)], [], 1.0, 2.0)
    model = residue_stats_model(vars, [con])

    result = PC.residue_presolve_constraint!(model, con, [2])
    stats = result.residue_stats

    @test result.changed
    @test result.tightened_to_equality
    @test con.lhs == 2.0
    @test con.rhs == 2.0
    @test stats.num_constraints_processed == 1
    @test stats.num_constraints_tightened == 1
    @test stats.num_constraints_eq_tightened == 1
    @test stats.num_lower_bounds_tightened == 1
    @test stats.num_upper_bounds_tightened == 0
    @test stats.num_moduli_evaluated == 1
    @test stats.num_useful_moduli == 1
    @test stats.avg_relative_bound_tightening ≈ 1.0
end

@testset "ResidueStats records residue-detected infeasibility" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 10.0))
    con = residue_stats_constraint([], [(2.0, 1)], 1.0, 1.0)
    model = residue_stats_model(vars, [con])

    result = PC.residue_presolve_constraint!(model, con, [2])
    stats = result.residue_stats

    @test result.changed
    @test result.processed
    @test result.infeasible
    @test model.infeasible
    @test con.lhs > con.rhs
    @test stats.infeasibility_detected
    @test stats.num_constraints_processed == 1
    @test stats.num_constraints_tightened == 1
    @test stats.num_constraints_eq_tightened == 0
    @test stats.num_lower_bounds_tightened == 1
    @test stats.num_upper_bounds_tightened == 1
    @test stats.num_moduli_evaluated == 1
    @test stats.num_useful_moduli == 1
    @test stats.avg_relative_bound_tightening == 0.0
end

@testset "ResidueStats records saturated non-useful moduli" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 5.0))
    con = residue_stats_constraint([], [(1.0, 1)], 0.0, 5.0)
    model = residue_stats_model(vars, [con])

    result = PC.residue_presolve_constraint!(model, con, [2])
    stats = result.residue_stats

    @test !result.changed
    @test result.processed
    @test !result.infeasible
    @test con.lhs == 0.0
    @test con.rhs == 5.0
    @test stats.num_constraints_processed == 1
    @test stats.num_constraints_tightened == 0
    @test stats.num_moduli_evaluated == 1
    @test stats.num_useful_moduli == 0
    @test stats.avg_relative_bound_tightening == 0.0
end

@testset "ResidueStats counts constraints uniquely with shared accumulator" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(1.0, 4.0))
    con = residue_stats_constraint([], [(2.0, 1)], 1.0, 9.0)
    model = residue_stats_model(vars, [con])
    accumulator = PC._ResidueStatsAccumulator()

    first = PC.residue_presolve_constraint!(model, con, [2], accumulator)
    second = PC.residue_presolve_constraint!(model, con, [2], accumulator)
    stats = accumulator.stats

    @test first.residue_stats === stats
    @test second.residue_stats === stats
    @test first.changed
    @test !second.changed
    @test stats.num_constraints_processed == 1
    @test stats.num_constraints_tightened == 1
    @test stats.num_lower_bounds_tightened == 1
    @test stats.num_upper_bounds_tightened == 1
    @test stats.num_moduli_evaluated == 2
    @test stats.num_useful_moduli == 1
    @test stats.avg_relative_bound_tightening ≈ 0.25
end

@testset "PresolveResult exposes accumulated ResidueStats" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 4.0),
        2 => PC.IntVar(0.0, 4.0),
    )
    con = residue_stats_constraint([(2.0, 1, 2)], [], 1.0, 9.0)
    model = residue_stats_model(vars, [con])

    result = QIPresolve.presolve!(
        model;
        residue_strategy = :powers_of_two,
        residue_threshold = 2,
        treewidth_threshold = 2,
    )
    stats = result.residue_stats

    @test result isa QIPresolve.PresolveResult
    @test result.parity_stats isa PC.ParityStats
    @test stats isa PC.ResidueStats
    @test result.model === model
    @test !model.infeasible
    @test model.cons[1].lhs == 2.0
    @test model.cons[1].rhs == 8.0
    @test stats.num_constraints_processed == 1
    @test stats.num_constraints_tightened == 1
    @test stats.num_lower_bounds_tightened == 1
    @test stats.num_upper_bounds_tightened == 1
    @test stats.num_moduli_evaluated == 1
    @test stats.num_useful_moduli == 1
    @test stats.avg_relative_bound_tightening ≈ 0.25
    @test stats.residue_presolve_time >= 0.0

    compatibility_result = QIPresolve.PresolveResult(
        residue_stats_model(Dict{PC.VarId, PC.IntVar}(), PC.Constraint[]),
        PC.ParityPostsolver(PC.VarId[]),
    )
    @test compatibility_result.parity_stats isa PC.ParityStats
    @test compatibility_result.residue_stats isa PC.ResidueStats

    compatibility_with_parity_stats = QIPresolve.PresolveResult(
        residue_stats_model(Dict{PC.VarId, PC.IntVar}(), PC.Constraint[]),
        PC.ParityPostsolver(PC.VarId[]),
        PC.ParityStats(),
    )
    @test compatibility_with_parity_stats.residue_stats isa PC.ResidueStats
end
