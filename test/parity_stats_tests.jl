using Test
import QIPresolve
import QIPresolve.PresolvingCore as PC

const ParityStatsQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const ParityStatsLinTerm = Tuple{Float64, PC.VarId}
const PARITY_STATS_NEXT_CON_ID = Ref(0)

parity_stats_next_con_id() = (PARITY_STATS_NEXT_CON_ID[] += 1)
parity_stats_empty_objective() = PC.QuadExpr(ParityStatsQuadTerm[], ParityStatsLinTerm[])

function parity_stats_model(vars, cons)
    return PC.QPModel(vars, cons, parity_stats_empty_objective(), :min)
end

function parity_stats_constraint(quad_terms, lin_terms, lhs::Float64, rhs::Float64)
    return PC.Constraint(
        parity_stats_next_con_id(),
        PC.QuadExpr(ParityStatsQuadTerm[quad_terms...], ParityStatsLinTerm[lin_terms...]),
        lhs,
        rhs,
    )
end

function parity_stats_linearized_fixed_model()
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 4.0),
        2 => PC.IntVar(0.0, 5.0),
        3 => PC.IntVar(0.0, 5.0),
        4 => PC.IntVar(0.0, 5.0),
    )
    fixed_parity_con = parity_stats_constraint(
        [(1.0, 1, 1)],
        [],
        0.0,
        0.0,
    )
    linearized_con = parity_stats_constraint(
        [(2.0, 1, 2)],
        [(1.0, 3), (1.0, 4)],
        0.0,
        10.0,
    )
    return parity_stats_model(vars, [fixed_parity_con, linearized_con])
end

function parity_stats_pattern_model()
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = parity_stats_constraint(
        [],
        [(1.0, 1), (1.0, 2)],
        4.0,
        4.0,
    )
    return parity_stats_model(vars, [con])
end

function parity_stats_single_distance_model()
    vars = Dict{PC.VarId, PC.IntVar}(i => PC.IntVar(-20.0, 20.0) for i in 1:4)
    con = parity_stats_constraint(
        [
            (1.0, 1, 1),
            (1.0, 2, 2),
            (-2.0, 1, 2),
            (1.0, 3, 3),
            (1.0, 4, 4),
            (-2.0, 3, 4),
        ],
        [],
        402.0,
        402.0,
    )
    return parity_stats_model(vars, [con])
end

@testset "ParityStats defaults and no-op parity presolve" begin
    defaults = PC.ParityStats()

    for field in fieldnames(PC.ParityStats)
        value = getfield(defaults, field)
        if value isa Integer
            @test value == 0
        end
    end
    @test defaults.parity_presolve_time == 0.0
    @test !defaults.infeasibility_detected
    @test defaults.infeasibility_source == :none

    @test isdefined(PC, :ParityStats)
    @test !isdefined(QIPresolve, :ParityStats)

    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(0.0, 5.0),
    )
    con = parity_stats_constraint(
        [],
        [(1.0, 1), (1.0, 2)],
        -Inf,
        3.0,
    )
    model = parity_stats_model(vars, [con])

    result = PC.parity_presolve!(model)
    stats = result.parity_stats

    @test stats isa PC.ParityStats
    @test !result.changed
    @test !model.infeasible
    @test stats.num_parity_presolve_rounds == 1
    @test stats.num_parity_presolve_phases == 1
    @test stats.num_parity_constraints_generated == 0
    @test stats.num_xor_constraints_generated == 0
    @test stats.num_xorand_constraints_generated == 0
    @test stats.num_parity_substitutions == 0
    @test stats.parity_presolve_time >= 0.0
    @test !stats.infeasibility_detected
    @test stats.infeasibility_source == :none
end

@testset "ParityStats records generated rows and canonical conjunction variables" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
    )
    con = parity_stats_constraint(
        [(2.0, 1, 2)],
        [],
        0.0,
        0.0,
    )
    model = parity_stats_model(vars, [con])
    acc = PC._ParityStatsAccumulator()
    parity_model = PC.build_parity_model(model, acc)
    stats = acc.stats

    @test length(parity_model.cons) == 2
    @test stats.num_parity_constraints_generated == 2
    @test stats.num_xor_constraints_generated == 1
    @test stats.num_xorand_constraints_generated == 1
    @test stats.num_parity_variables == 2
    @test stats.num_conjunction_variables == 1
    @test stats.num_total_boolean_variables == 3

    permuted_qe = PC.QuadExpr(
        ParityStatsQuadTerm[(0.0, 1, 99), (2.0, 1, 2)],
        ParityStatsLinTerm[],
    )
    @test PC.remove_var!(permuted_qe, 99)
    @test PC.get_quad_coeff(permuted_qe, 1, 2) == 2.0

    permuted_model = parity_stats_model(
        deepcopy(vars),
        [PC.Constraint(parity_stats_next_con_id(), permuted_qe, 0.0, 0.0)],
    )
    permuted_acc = PC._ParityStatsAccumulator()
    PC.build_parity_model(permuted_model, permuted_acc)

    @test permuted_acc.stats.num_xorand_constraints_generated == 1
    @test permuted_acc.stats.num_conjunction_variables == 1
    @test permuted_acc.stats.num_total_boolean_variables == 3
end

@testset "ParityStats records propagation edges and fixing sources" begin
    edge_acc = PC._ParityStatsAccumulator()
    edge_manager = PC.PropagationManager([1, 2])

    PC.add_implication!(edge_manager, PC.VarLit(1, false), PC.VarLit(2, false), edge_acc)
    PC.add_implication!(edge_manager, PC.VarLit(1, false), PC.VarLit(2, false), edge_acc)
    @test edge_acc.stats.num_implication_edges_generated == 2

    elimination_acc = PC._ParityStatsAccumulator()
    elimination_manager = PC.PropagationManager([1])
    PC.fix_var!(elimination_manager, 1, true, elimination_acc, :elimination)
    @test elimination_acc.stats.num_parity_fixed_elimination == 1
    @test elimination_acc.stats.num_implication_edges_generated == 1

    propagation_acc = PC._ParityStatsAccumulator()
    propagation_manager = PC.PropagationManager([1])
    @test PC.set_lit_label!(propagation_manager, 1, PC.TRUE, propagation_acc, :propagation)
    @test PC.set_lit_label!(propagation_manager, 2, PC.FALSE, propagation_acc, :propagation)
    @test propagation_acc.stats.num_parity_fixed_propagation == 1

    path_acc = PC._ParityStatsAccumulator()
    path_manager = PC.PropagationManager([1])
    PC.add_edge!(path_manager, PC.VarLit(1, true), PC.VarLit(1, false), path_acc)
    PC.update!(path_manager, path_acc)
    @test path_acc.stats.num_parity_fixed_path == 1
    @test path_acc.stats.num_implication_edges_generated == 1
end

@testset "ParityStats records fixed parity substitution effects" begin
    model = parity_stats_linearized_fixed_model()

    result = PC.parity_presolve!(model)
    stats = result.parity_stats

    @test result.fixed_parities == 2
    @test !model.infeasible
    @test length(model.cons) == 1
    @test stats.num_parity_constraints_generated == 4
    @test stats.num_xor_constraints_generated == 4
    @test stats.num_parity_fixed_elimination == 2
    @test stats.num_fixed_parity_substitutions == 2
    @test stats.num_pattern_substitutions == 0
    @test stats.num_parity_substitutions == 2
    @test stats.num_variables_fixed_after_parity_substitution == 1
    @test stats.num_constraints_modified_by_parity == 1
    @test stats.num_quadratic_terms_removed_by_parity == 2
    @test stats.num_constraints_linearized_by_parity == 1
    @test stats.num_constraints_removed_by_parity == 1
    @test stats.num_parity_presolve_rounds == 1
    @test stats.num_parity_presolve_phases == 3
end

@testset "ParityStats records pattern substitutions and split xor-and rows" begin
    pattern_model = parity_stats_pattern_model()

    pattern_result = PC.parity_presolve!(pattern_model)
    pattern_stats = pattern_result.parity_stats

    @test pattern_result.pattern_rewritten_vars == 2
    @test pattern_stats.num_pattern_relations == 1
    @test pattern_stats.num_pattern_substitutions == 2
    @test pattern_stats.num_parity_substitutions == 2
    @test pattern_stats.num_constraints_modified_by_parity == 1
    @test pattern_stats.num_quadratic_terms_removed_by_parity == 0
    @test pattern_stats.num_constraints_linearized_by_parity == 0
    @test pattern_stats.num_constraints_removed_by_parity == 0

    single_distance_model = parity_stats_single_distance_model()
    build_acc = PC._ParityStatsAccumulator()
    parity_model = PC.build_parity_model(single_distance_model, build_acc)

    @test length(parity_model.cons) == 2
    @test build_acc.stats.num_parity_constraints_generated == 2
    @test build_acc.stats.num_xor_constraints_generated == 1
    @test build_acc.stats.num_xorand_constraints_generated == 1
    @test build_acc.stats.num_parity_variables == 4
    @test build_acc.stats.num_conjunction_variables == 4
    @test build_acc.stats.num_total_boolean_variables == 8

    split_result = PC.parity_presolve!(single_distance_model)
    split_stats = split_result.parity_stats

    @test split_result.pattern_rewritten_vars == 4
    @test split_stats.num_constraints_split == 1
    @test split_stats.num_parity_constraints_generated == 5
    @test split_stats.num_xor_constraints_generated == 4
    @test split_stats.num_xorand_constraints_generated == 1
    @test split_stats.num_pattern_relations == 2
    @test split_stats.num_pattern_substitutions == 4
    @test split_stats.num_parity_substitutions == 4
end

@testset "ParityStats records parity infeasibility sources" begin
    elimination_acc = PC._ParityStatsAccumulator()
    elimination_model = PC.ParityModel(
        Dict{PC.VarId, Int}(1 => 1, 2 => 2),
        [1, 2],
        PC.XorConstraint[
            PC.XorConstraint(BitVector([1, 1]), false),
            PC.XorConstraint(BitVector([1, 1]), true),
        ],
    )

    PC.gauss_jordan_xor!(elimination_model, elimination_acc)
    @test elimination_model.infeasible
    @test elimination_acc.stats.infeasibility_detected
    @test elimination_acc.stats.infeasibility_source == :elimination

    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
        3 => PC.IntVar(0.0, 1.0),
    )
    con1 = parity_stats_constraint([], [(1.0, 1), (1.0, 2)], 0.0, 0.0)
    con2 = parity_stats_constraint([], [(1.0, 2), (1.0, 3)], 0.0, 0.0)
    con3 = parity_stats_constraint([], [(1.0, 1), (1.0, 3)], 1.0, 1.0)
    propagation_model = parity_stats_model(vars, [con1, con2, con3])

    result = PC.parity_presolve!(propagation_model)
    stats = result.parity_stats

    @test propagation_model.infeasible
    @test result.infeasible
    @test stats.infeasibility_detected
    @test stats.infeasibility_source == :propagation
    @test stats.num_parity_constraints_generated == 6
    @test stats.num_implication_edges_generated == 16

    entry_model = parity_stats_model(Dict{PC.VarId, PC.IntVar}(), PC.Constraint[])
    entry_model.infeasible = true
    entry_result = PC.parity_presolve!(entry_model)

    @test entry_result.infeasible
    @test entry_result.parity_stats.num_parity_presolve_rounds == 0
    @test entry_result.parity_stats.num_parity_presolve_phases == 0
    @test !entry_result.parity_stats.infeasibility_detected
    @test entry_result.parity_stats.infeasibility_source == :none
end

@testset "PresolveResult exposes accumulated ParityStats" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
    )
    con = parity_stats_constraint(
        [(2.0, 1, 2)],
        [],
        1.0,
        2.0,
    )
    model = parity_stats_model(vars, [con])

    result = QIPresolve.presolve!(
        model;
        residue_strategy = :powers_of_two,
        residue_threshold = 2,
    )
    stats = result.parity_stats

    @test result isa QIPresolve.PresolveResult
    @test stats isa PC.ParityStats
    @test stats.num_parity_presolve_rounds == 2
    @test stats.num_parity_presolve_phases == 3
    @test stats.num_parity_constraints_generated == 2
    @test stats.num_xor_constraints_generated == 1
    @test stats.num_xorand_constraints_generated == 1
    @test stats.num_fixed_parity_substitutions == 2
    @test stats.num_parity_substitutions == 2
    @test stats.parity_presolve_time >= 0.0

    compatibility_result = QIPresolve.PresolveResult(
        parity_stats_model(Dict{PC.VarId, PC.IntVar}(), PC.Constraint[]),
        PC.ParityPostsolver(PC.VarId[]),
    )

    @test compatibility_result.parity_stats isa PC.ParityStats
    @test compatibility_result.parity_stats.num_parity_presolve_rounds == 0
end
