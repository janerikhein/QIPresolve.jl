using Test
import QIPresolve.PresolvingCore as PC

const ParityQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const ParityLinTerm = Tuple{Float64, PC.VarId}
const PARITY_NEXT_CON_ID = Ref(0)

parity_next_con_id() = (PARITY_NEXT_CON_ID[] += 1)

function parity_empty_objective()
    return PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[])
end

function parity_edge_matrix(n::Int, edges::Vector{Tuple{Int, Int}})
    mat = falses(n, n)
    for (i, j) in edges
        mat[i, j] = true
        mat[j, i] = true
    end
    return mat
end

function parity_single_distance_model()
    vars = Dict{PC.VarId, PC.IntVar}(i => PC.IntVar(-20.0, 20.0) for i in 1:4)
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(
            ParityQuadTerm[
                (1.0, 1, 1),
                (1.0, 2, 2),
                (-2.0, 1, 2),
                (1.0, 3, 3),
                (1.0, 4, 4),
                (-2.0, 3, 4),
            ],
            ParityLinTerm[],
        ),
        402.0,
        402.0,
    )
    return PC.QPModel(vars, [con], parity_empty_objective(), :min)
end

function parity_eval_equivalent(model0, model, var_id, scale, offset)
    con0 = model0.cons[1]
    con = model.cons[1]
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs

    @test isapprox(shift_lhs, shift_rhs; atol = 1.0e-12)

    for _ in 1:5
        x = randn(2)
        x_sub = copy(x)
        x_sub[var_id] = scale * x[var_id] + offset

        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1.0e-8)

        obj_before = PC.eval_full(model0.obj_expr, x_sub)
        obj_after = PC.eval_full(model.obj_expr, x)
        @test isapprox(obj_after, obj_before; atol = 1.0e-8)
    end
end

function parity_prepare_manager(lit1::PC.VarLit, lit2::PC.VarLit)
    manager = PC.PropagationManager([1, 2])
    PC.add_equivalence!(manager, lit1, lit2)
    PC.update_sccs!(manager)
    while PC.pop_substitution!(manager) !== nothing
    end
    return manager
end

@testset "count_builder_occurrences! ignores toggled-off builder entries" begin
    counts = Dict{PC.VarId, Int}()
    builder = PC.XorConstraintBuilder()

    PC.add_par!(builder, 7)
    PC.add_par!(builder, 7)
    PC.add_par!(builder, 3)

    PC.add_conj!(builder, 2, 5)
    PC.add_conj!(builder, 5, 2)
    PC.add_conj!(builder, 3, 5)

    PC.count_builder_occurrences!(counts, builder)

    @test counts == Dict(3 => 2, 5 => 1)
    @test !haskey(counts, 2)
    @test !haskey(counts, 7)
end

@testset "build_parity_model orders active vars by occurrence count and preserves row order" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(-1.0, 2.0),
        2 => PC.IntVar(-1.0, 2.0),
        3 => PC.IntVar(0.0, 1.0),
        4 => PC.IntVar(0.0, 1.0),
    )

    con1 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 2)]),
        0.0,
        0.0,
    )
    con2 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 3)]),
        0.0,
        0.0,
    )
    con3 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 2)], ParityLinTerm[]),
        0.0,
        0.0,
    )
    con4 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 2, 3)], ParityLinTerm[]),
        0.0,
        0.0,
    )

    model = PC.QPModel(vars, [con1, con2, con3, con4], parity_empty_objective(), :min)
    parity_model = PC.build_parity_model(model)
    mod2_model = PC.build_parity_model(model; parity_strategy = :mod2_basic)
    mod4_basic_model = PC.build_parity_model(model; parity_strategy = :mod4_basic)
    full_model = PC.build_parity_model(model; parity_strategy = :full)

    @test parity_model.pos_to_var_id == [2, 1, 3]
    @test parity_model.var_id_to_pos == Dict(2 => 1, 1 => 2, 3 => 3)
    @test !haskey(parity_model.var_id_to_pos, 4)
    @test length(parity_model.cons) == 6
    @test length(mod2_model.cons) == 4
    @test all(con -> con.meta.is_pure_xor, mod2_model.cons)
    @test length(mod4_basic_model.cons) == 6
    @test length(full_model.cons) == 6
    @test_throws ArgumentError PC.build_parity_model(model; parity_strategy = :invalid)

    @test parity_model.cons[1].meta.is_pure_xor
    @test parity_model.cons[1].par == BitVector([1, 0, 0])
    @test !parity_model.cons[1].rhs

    @test parity_model.cons[2].meta.is_pure_xor
    @test parity_model.cons[2].par == BitVector([0, 1, 1])
    @test !parity_model.cons[2].rhs

    @test parity_model.cons[3].meta.is_pure_xor
    @test parity_model.cons[3].par == falses(3)
    @test parity_model.cons[3].conj === nothing

    @test !parity_model.cons[4].meta.is_pure_xor
    @test parity_model.cons[4].par == falses(3)
    @test parity_model.cons[4].conj == parity_edge_matrix(3, [(1, 2)])

    @test parity_model.cons[5].meta.is_pure_xor
    @test parity_model.cons[5].par == falses(3)
    @test parity_model.cons[5].conj === nothing

    @test !parity_model.cons[6].meta.is_pure_xor
    @test parity_model.cons[6].par == falses(3)
    @test parity_model.cons[6].conj == parity_edge_matrix(3, [(1, 3)])
end


@testset "parity XOR branch expands conjunction terms using pivot substitutions" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
        3 => PC.IntVar(0.0, 1.0),
        4 => PC.IntVar(0.0, 1.0),
    )

    xor_con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2), (1.0, 3)]),
        1.0,
        1.0,
    )
    mixed_con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 4)], ParityLinTerm[(2.0, 4)]),
        0.0,
        0.0,
    )

    model = PC.QPModel(vars, [xor_con, mixed_con], parity_empty_objective(), :min)
    parity_model = PC.build_parity_model(model)
    propagator = PC.PropagationManager(parity_model.pos_to_var_id)

    PC.propagate!(parity_model, propagator)
    PC.gauss_jordan_xor!(parity_model)
    PC.propagate!(parity_model, propagator)
    PC.substitute_pivots_in_conjunctive_terms!(parity_model)

    @test parity_model.pos_to_var_id == [1, 2, 3, 4]
    @test length(parity_model.cons) == 3
    @test parity_model.pivots[1] == (1, nothing)
    @test parity_model.cons[1].par == BitVector([1, 1, 1, 0])
    @test parity_model.cons[1].rhs
    @test parity_model.cons[2].par == falses(4)
    @test parity_model.cons[2].conj == parity_edge_matrix(4, [(2, 3)])
    @test !parity_model.cons[2].rhs
    @test parity_model.cons[3].par == falses(4)
    @test parity_model.cons[3].conj == parity_edge_matrix(4, [(2, 4), (3, 4)])
    @test !parity_model.cons[3].rhs
end

@testset "parity XOR branch eliminates stored parity pivots from xor-and rows" begin
    pos_to_var_id = [1, 2, 3, 4]
    var_id_to_pos = Dict(i => i for i in pos_to_var_id)
    xor_pivot = PC.XorConstraint(BitVector([1, 1, 0, 0]), false)
    xor_and_row = PC.XorConstraint(
        BitVector([1, 0, 0, 0]),
        parity_edge_matrix(4, [(3, 4)]),
        true,
    )
    parity_model = PC.ParityModel(var_id_to_pos, pos_to_var_id, [xor_pivot, xor_and_row])
    parity_model.pivots[1] = (1, nothing)

    PC.substitute_pivots_in_conjunctive_terms!(parity_model)
    PC.substitute_parity_pivots!(parity_model)

    @test parity_model.cons[2].par == BitVector([0, 1, 0, 0])
    @test parity_model.cons[2].conj == parity_edge_matrix(4, [(3, 4)])
    @test parity_model.cons[2].rhs
end

@testset "fix_parities! rewrites even and odd parity-fixed variables" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(-2.0, 5.0),
        2 => PC.IntVar(0.0, 3.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(
            ParityQuadTerm[(1.5, 1, 1), (-0.75, 1, 2)],
            ParityLinTerm[(2.0, 1), (0.5, 2)];
            constant = -1.0,
        ),
        -4.0,
        6.0,
    )
    obj = PC.QuadExpr(
        ParityQuadTerm[(1.25, 1, 2)],
        ParityLinTerm[(-3.0, 1), (1.0, 2)];
        constant = 2.0,
    )

    even_model0 = PC.QPModel(deepcopy(vars), [deepcopy(con)], deepcopy(obj), :min)
    even_model = deepcopy(even_model0)
    even_manager = PC.PropagationManager([1])
    PC.fix_var!(even_manager, 1, false)

    @test PC.fix_parities!(even_model, even_manager) == 1
    @test even_model.vars[1] == PC.IntVar(-1.0, 2.0)
    parity_eval_equivalent(even_model0, even_model, 1, 2.0, 0.0)

    odd_model0 = PC.QPModel(deepcopy(vars), [deepcopy(con)], deepcopy(obj), :min)
    odd_model = deepcopy(odd_model0)
    odd_manager = PC.PropagationManager([1])
    PC.fix_var!(odd_manager, 1, true)

    @test PC.fix_parities!(odd_model, odd_manager) == 1
    @test odd_model.vars[1] == PC.IntVar(-1.0, 2.0)
    parity_eval_equivalent(odd_model0, odd_model, 1, 2.0, 1.0)
end

@testset "fix_parities! is a no-op when no fixed parity is available" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 4.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1)]),
        0.0,
        2.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    model0 = deepcopy(model)
    manager = PC.PropagationManager([1])

    @test PC.fix_parities!(model, manager) == 0
    @test model.vars == model0.vars
    @test model.cons[1].lhs == model0.cons[1].lhs
    @test model.cons[1].rhs == model0.cons[1].rhs
    @test model.cons[1].qe.constant == model0.cons[1].qe.constant
    @test collect(PC.vars(model.cons[1])) == collect(PC.vars(model0.cons[1]))
end

@testset "fix_parities! ignores parity fixings for absent variables" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 4.0))
    model = PC.QPModel(vars, PC.Constraint[], parity_empty_objective(), :min)
    model0 = deepcopy(model)
    manager = PC.PropagationManager([99])
    PC.fix_var!(manager, 99, true)

    @test PC.fix_parities!(model, manager) == 0
    @test model.vars == model0.vars
    @test isempty(model.cons)
end

@testset "fix_parities! keeps singleton transformed variables in the model" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1)]),
        1.0,
        1.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    manager = PC.PropagationManager([1])
    PC.fix_var!(manager, 1, true)

    @test PC.fix_parities!(model, manager) == 1
    @test haskey(model.vars, 1)
    @test model.vars[1] == PC.IntVar(0.0, 0.0)
end

@testset "fix_parity_patterns! introduces one binary variable for a positive SCC" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        0.0,
        10.0,
    )
    obj = PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (-2.0, 2)])
    model = PC.QPModel(vars, [con], obj, :min)
    manager = parity_prepare_manager(PC.VarLit(1, false), PC.VarLit(2, false))

    @test PC.fix_parity_patterns!(model, manager) == 2
    @test haskey(model.vars, 3)
    @test PC.is_binary(model.vars[3])
    @test model.vars[1] == PC.IntVar(0.0, 2.0)
    @test model.vars[2] == PC.IntVar(0.0, 3.0)
    @test length(model.cons) == 5

    pattern_cons = model.cons[2:end]
    @test pattern_cons[1].lhs == 0.0
    @test pattern_cons[1].rhs == 4.0
    @test PC.get_lin_coeff(pattern_cons[1].qe, 1) == 1.0
    @test PC.get_lin_coeff(pattern_cons[1].qe, 3) == 2.0
    @test pattern_cons[2].lhs == -2.0
    @test pattern_cons[2].rhs == 2.0
    @test PC.get_lin_coeff(pattern_cons[2].qe, 1) == 1.0
    @test PC.get_lin_coeff(pattern_cons[2].qe, 3) == -2.0
    @test pattern_cons[3].lhs == 1.0
    @test pattern_cons[3].rhs == 5.0
    @test PC.get_lin_coeff(pattern_cons[3].qe, 2) == 1.0
    @test PC.get_lin_coeff(pattern_cons[3].qe, 3) == 3.0
    @test pattern_cons[4].lhs == -3.0
    @test pattern_cons[4].rhs == 3.0
    @test PC.get_lin_coeff(pattern_cons[4].qe, 2) == 1.0
    @test PC.get_lin_coeff(pattern_cons[4].qe, 3) == -3.0

    @test manager.lit_to_pos[PC.VarLit(3, false)] == manager.scc_pos_to_rep_pos[manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(3, false)]]]
    @test manager.lit_to_pos[PC.VarLit(3, true)] == manager.scc_pos_to_rep_pos[manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(3, true)]]]
    @test PC.fixed_value(manager, 3) === nothing
end

@testset "fix_parity_patterns! handles negated literals in the positive SCC" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        0.0,
        10.0,
    )
    obj = PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (-2.0, 2)])
    model = PC.QPModel(vars, [con], obj, :min)
    model0 = deepcopy(model)
    manager = parity_prepare_manager(PC.VarLit(1, false), PC.VarLit(2, true))

    @test PC.fix_parity_patterns!(model, manager) == 2
    @test haskey(model.vars, 3)
    @test PC.is_binary(model.vars[3])
    @test model.vars[1] == PC.IntVar(0.0, 2.0)
    @test model.vars[2] == PC.IntVar(0.0, 3.0)
    @test length(model.cons) == 5

    pattern_cons = model.cons[2:end]
    @test pattern_cons[3].lhs == 0.0
    @test pattern_cons[3].rhs == 6.0
    @test PC.get_lin_coeff(pattern_cons[3].qe, 2) == 1.0
    @test PC.get_lin_coeff(pattern_cons[3].qe, 3) == 3.0
    @test pattern_cons[4].lhs == -2.0
    @test pattern_cons[4].rhs == 2.0
    @test PC.get_lin_coeff(pattern_cons[4].qe, 2) == 1.0
    @test PC.get_lin_coeff(pattern_cons[4].qe, 3) == -3.0

    for z in (0.0, 1.0)
        for y1 in 0.0:2.0
            for y2 in 0.0:3.0
                x = zeros(3)
                x[1] = y1
                x[2] = y2
                x[3] = z

                x0 = zeros(2)
                x0[1] = 2.0 * y1 + z
                x0[2] = 2.0 * y2 - z + 1.0

                val_before = PC.eval_full(model0.cons[1].qe, x0)
                val_after = PC.eval_full(model.cons[1].qe, x)
                shift = model0.cons[1].lhs - model.cons[1].lhs
                @test isapprox(val_after, val_before - shift; atol = 1.0e-8)

                obj_before = PC.eval_full(model0.obj_expr, x0)
                obj_after = PC.eval_full(model.obj_expr, x)
                @test isapprox(obj_after, obj_before; atol = 1.0e-8)
            end
        end
    end
end

@testset "fix_parity_patterns! is a no-op without an eligible SCC" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 5.0))
    model = PC.QPModel(vars, PC.Constraint[], parity_empty_objective(), :min)
    model0 = deepcopy(model)
    manager = PC.PropagationManager([1])

    @test PC.fix_parity_patterns!(model, manager) == 0
    @test model.vars == model0.vars
    @test isempty(model.cons)
end

@testset "parity_presolve! removes fixed parity singletons during normalization" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1)]),
        1.0,
        1.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)

    stats = PC.parity_presolve!(model)

    @test stats.changed
    @test stats.domains_changed
    @test isempty(stats.coefficient_changed_constraint_ids)
    @test stats.fixed_parities == 0
    @test stats.pattern_rewritten_vars == 0
    @test !stats.infeasible
    @test !model.infeasible
    @test !haskey(model.vars, 1)
    @test isempty(model.cons)
end

@testset "parity_presolve! scales constraints by gcd after parity rewrites" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 4.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(1.0, 1, 1)], ParityLinTerm[]),
        0.0,
        0.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)

    PC.parity_presolve!(model)

    @test !model.infeasible
    @test isempty(model.vars)
    @test isempty(model.cons)
end

@testset "parity_presolve! handles negative singleton constraints produced by rewrites" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 4.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(1.0, 1, 1)], ParityLinTerm[]),
        0.0,
        0.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)

    PC.parity_presolve!(model)

    @test !model.infeasible
end

@testset "parity_presolve_phase! returns whether a phase changed the QP model" begin
    changed_vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    changed_con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        4.0,
        4.0,
    )
    changed_model = PC.QPModel(changed_vars, [changed_con], parity_empty_objective(), :min)
    propagator = PC.PropagationManager(PC.VarId[])

    changed_stats = PC.parity_presolve_phase!(changed_model, propagator)
    @test changed_stats == (changed = true, fixed_parities = 0, pattern_rewritten_vars = 2)
    @test !changed_model.infeasible
    @test haskey(changed_model.vars, 3)
    @test PC.is_binary(changed_model.vars[3])
    @test changed_model.vars[1] == PC.IntVar(0.0, 2.0)
    @test changed_model.vars[2] == PC.IntVar(0.0, 3.0)

    stable_vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(0.0, 5.0),
    )
    stable_con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        -Inf,
        3.0,
    )
    stable_model = PC.QPModel(stable_vars, [stable_con], parity_empty_objective(), :min)

    stable_stats = PC.parity_presolve_phase!(stable_model, propagator)
    @test stable_stats == (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)
    @test !stable_model.infeasible
    @test stable_model.vars == stable_vars
    @test length(stable_model.cons) == 1
end

@testset "parity_presolve! introduces parity patterns and reaches a fixed point" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        4.0,
        4.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)

    PC.parity_presolve!(model)

    @test !model.infeasible
    @test haskey(model.vars, 3)
    @test PC.is_binary(model.vars[3])
    @test model.vars[1] == PC.IntVar(0.0, 2.0)
    @test model.vars[2] == PC.IntVar(0.0, 3.0)
    @test length(model.cons) == 5

    parity_model = PC.build_parity_model(model)
    propagator = PC.PropagationManager(parity_model.pos_to_var_id)
    PC.propagate!(parity_model, propagator)
    @test !parity_model.infeasible
end

@testset "parity_presolve! splits the single-distance xor-and row and halves all domains" begin
    extraction_model = parity_single_distance_model()
    parity_model = PC.build_parity_model(extraction_model)

    @test parity_model.pos_to_var_id == [1, 2, 3, 4]
    @test parity_model.var_id_to_pos == Dict(i => i for i in 1:4)
    @test length(parity_model.cons) == 2

    @test parity_model.cons[1].meta.is_pure_xor
    @test parity_model.cons[1].par == BitVector([1, 1, 1, 1])
    @test !parity_model.cons[1].rhs

    @test !parity_model.cons[2].meta.is_pure_xor
    @test parity_model.cons[2].par == falses(4)
    @test parity_model.cons[2].conj == parity_edge_matrix(4, [(1, 3), (1, 4), (2, 3), (2, 4)])
    @test parity_model.cons[2].rhs

    phase_model = parity_single_distance_model()
    propagator = PC.PropagationManager(PC.VarId[])
    phase_stats = PC.parity_presolve_phase!(phase_model, propagator)

    @test phase_stats == (changed = true, fixed_parities = 0, pattern_rewritten_vars = 4)
    @test !phase_model.infeasible
    @test phase_model.vars[1] == PC.IntVar(-10.0, 10.0)
    @test phase_model.vars[2] == PC.IntVar(-10.0, 10.0)
    @test phase_model.vars[3] == PC.IntVar(-10.0, 10.0)
    @test phase_model.vars[4] == PC.IntVar(-10.0, 10.0)
    @test haskey(phase_model.vars, 5)
    @test haskey(phase_model.vars, 6)
    @test PC.is_binary(phase_model.vars[5])
    @test PC.is_binary(phase_model.vars[6])

    full_model = parity_single_distance_model()
    PC.parity_presolve!(full_model)

    @test !full_model.infeasible
    @test sort!(collect(keys(full_model.vars))) == [1, 2, 3, 4, 5, 6]
    @test full_model.vars[1] == PC.IntVar(-10.0, 10.0)
    @test full_model.vars[2] == PC.IntVar(-10.0, 10.0)
    @test full_model.vars[3] == PC.IntVar(-10.0, 10.0)
    @test full_model.vars[4] == PC.IntVar(-10.0, 10.0)
    @test PC.is_binary(full_model.vars[5])
    @test PC.is_binary(full_model.vars[6])
    @test length(full_model.cons) == 9
end

@testset "parity_presolve_phase! reuses one propagator across multiple phases" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        4.0,
        4.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    propagator = PC.PropagationManager(PC.VarId[])

    phase1_stats = PC.parity_presolve_phase!(model, propagator)
    @test phase1_stats == (changed = true, fixed_parities = 0, pattern_rewritten_vars = 2)
    @test haskey(model.vars, 3)
    @test haskey(propagator.lit_to_pos, PC.VarLit(3, false))
    @test haskey(propagator.lit_to_pos, PC.VarLit(3, true))
    first_manager_size = length(propagator.pos_to_lit)

    phase2_stats = PC.parity_presolve_phase!(model, propagator)
    @test phase2_stats == (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)
    @test !model.infeasible
    @test length(propagator.pos_to_lit) == first_manager_size
    @test haskey(propagator.lit_to_pos, PC.VarLit(3, false))
end

@testset "parity_presolve_phase! reports fixed parity counts" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 4.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(1.0, 1, 1)], ParityLinTerm[]),
        0.0,
        0.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    propagator = PC.PropagationManager(PC.VarId[])

    stats = PC.parity_presolve_phase!(model, propagator)

    @test stats == (changed = true, fixed_parities = 1, pattern_rewritten_vars = 0)
    @test !model.infeasible
    @test model.vars[1] == PC.IntVar(0.0, 2.0)
end

@testset "parity_presolve! marks models infeasible when parity propagation contradicts" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0))
    con1 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1)]),
        0.0,
        0.0,
    )
    con2 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1)]),
        1.0,
        1.0,
    )
    model = PC.QPModel(vars, [con1, con2], parity_empty_objective(), :min)

    stats = PC.parity_presolve!(model)

    @test model.infeasible
    @test stats.infeasible
end

@testset "parity_presolve! reports initially infeasible models" begin
    model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(),
        PC.Constraint[],
        parity_empty_objective(),
        :min,
    )
    model.infeasible = true

    stats = PC.parity_presolve!(model)

    @test !stats.changed
    @test !stats.domains_changed
    @test isempty(stats.coefficient_changed_constraint_ids)
    @test stats.fixed_parities == 0
    @test stats.pattern_rewritten_vars == 0
    @test stats.infeasible
end

@testset "parity_presolve! is a no-op without parity structure" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(0.0, 5.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        -Inf,
        3.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    model0 = deepcopy(model)

    PC.parity_presolve!(model)

    @test !model.infeasible
    @test model.vars == model0.vars
    @test length(model.cons) == 1
    @test model.cons[1].lhs == model0.cons[1].lhs
    @test model.cons[1].rhs == model0.cons[1].rhs
    @test collect(PC.vars(model.cons[1])) == collect(PC.vars(model0.cons[1]))
end

@testset "parity_presolve! aggregates parallel constraints at entry" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 2.0),
        2 => PC.IntVar(0.0, 2.0),
    )
    con1 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 2)], ParityLinTerm[]),
        0.0,
        8.0,
    )
    con2 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 2)], ParityLinTerm[]),
        2.0,
        6.0,
    )
    model = PC.QPModel(vars, [con1, con2], parity_empty_objective(), :min)

    stats = PC.parity_presolve!(model)

    @test stats.changed
    @test !model.infeasible
    @test length(model.cons) == 1
    @test model.cons[1].id == con1.id
    @test model.cons[1].lhs == 2.0
    @test model.cons[1].rhs == 6.0
end

@testset "parity_presolve_phase! skips parallel aggregation during phase normalization" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 2.0),
        2 => PC.IntVar(0.0, 2.0),
    )
    con1 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 2)], ParityLinTerm[]),
        0.0,
        8.0,
    )
    con2 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 2)], ParityLinTerm[]),
        2.0,
        6.0,
    )
    model = PC.QPModel(vars, [con1, con2], parity_empty_objective(), :min)
    propagator = PC.PropagationManager(PC.VarId[])

    stats = PC.parity_presolve_phase!(model, propagator)

    @test !stats.changed
    @test !model.infeasible
    @test length(model.cons) == 2
end

@testset "postsolve returns original variables for untouched mappings" begin
    postsolver = PC.ParityPostsolver([2, 1])
    reduced_solution = Dict{PC.VarId, Float64}(1 => 3.0, 2 => -1.0, 9 => 7.0)

    original_solution = PC.postsolve(postsolver, reduced_solution)

    @test original_solution == Dict{PC.VarId, Float64}(1 => 3.0, 2 => -1.0)
    @test original_solution !== reduced_solution
end

@testset "postsolve adds normalization offsets for shifted binaries" begin
    postsolver = PC.ParityPostsolver([1])
    PC.add_reconstruction_offset!(postsolver, 1, 3.0)

    original_solution = PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 1.0))

    @test original_solution == Dict{PC.VarId, Float64}(1 => 4.0)
end

@testset "postsolve rebuilds multiple fixed low-order bits in LSB order" begin
    postsolver = PC.ParityPostsolver([1])
    PC.append_fixed_bit!(postsolver, 1, true)
    PC.append_fixed_bit!(postsolver, 1, false)

    original_solution = PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 2.0))

    @test original_solution == Dict{PC.VarId, Float64}(1 => 9.0)
end

@testset "postsolve evaluates live binary bit references" begin
    postsolver = PC.ParityPostsolver([1, 2])
    PC.append_binary_bit!(postsolver, 1, 7)
    PC.append_binary_bit!(postsolver, 2, 7; negated = true)

    original_solution = PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 2.0, 2 => 2.0, 7 => 1.0))

    @test original_solution == Dict{PC.VarId, Float64}(1 => 5.0, 2 => 4.0)
end

@testset "fix_parities! and fix_vars! keep removed transformed values in the postsolver" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(1.0, 1.0))
    model = PC.QPModel(vars, PC.Constraint[], parity_empty_objective(), :min)
    manager = PC.PropagationManager([1])
    postsolver = PC.ParityPostsolver(keys(vars))
    PC.fix_var!(manager, 1, true)

    @test PC.fix_parities!(model, manager, postsolver) == 1
    PC.fix_vars!(model, postsolver)

    var_data = postsolver.var_data[1]
    @test !haskey(model.vars, 1)
    @test var_data.current_var_id === nothing
    @test var_data.fixed_high_order == 0.0
    @test length(var_data.bits) == 1
    @test var_data.bits[1].kind == PC.FIXED1
    @test var_data.bits[1].ref_var_id === nothing
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}()) == Dict{PC.VarId, Float64}(1 => 1.0)
end

@testset "fix_parity_patterns! postsolve resolves helper vars fixed later" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        0.0,
        10.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    manager = parity_prepare_manager(PC.VarLit(1, false), PC.VarLit(2, true))
    postsolver = PC.ParityPostsolver(keys(vars))

    @test PC.fix_parity_patterns!(model, manager, postsolver) == 2
    @test length(postsolver.var_data[1].bits) == 1
    @test postsolver.var_data[1].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[1].bits[1].ref_var_id == 3
    @test length(postsolver.var_data[2].bits) == 1
    @test postsolver.var_data[2].bits[1].kind == PC.NEGATED_BINVAR
    @test postsolver.var_data[2].bits[1].ref_var_id == 3

    PC.set_var_bounds!(model, 3, 1.0, 1.0)
    PC.fix_vars!(model, postsolver)

    @test !haskey(model.vars, 3)
    @test postsolver.var_data[3].current_var_id === nothing
    @test postsolver.var_data[3].fixed_high_order == 1.0
    @test isempty(postsolver.var_data[3].bits)
    @test postsolver.var_data[1].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[1].bits[1].ref_var_id == 3
    @test postsolver.var_data[2].bits[1].kind == PC.NEGATED_BINVAR
    @test postsolver.var_data[2].bits[1].ref_var_id == 3
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 2.0, 2 => 3.0)) ==
          Dict{PC.VarId, Float64}(1 => 5.0, 2 => 6.0)
end

@testset "fix_parity_patterns! postsolve resolves helper vars parity-fixed and removed later" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        0.0,
        10.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    manager = parity_prepare_manager(PC.VarLit(1, false), PC.VarLit(2, true))
    postsolver = PC.ParityPostsolver(keys(vars))

    @test PC.fix_parity_patterns!(model, manager, postsolver) == 2

    helper_manager = PC.PropagationManager([3])
    PC.fix_var!(helper_manager, 3, true)
    @test PC.fix_parities!(model, helper_manager, postsolver) == 1
    PC.fix_vars!(model, postsolver)

    @test !haskey(model.vars, 3)
    @test haskey(postsolver.var_data, 3)
    @test postsolver.var_data[3].current_var_id === nothing
    @test postsolver.var_data[3].fixed_high_order == 0.0
    @test length(postsolver.var_data[3].bits) == 1
    @test postsolver.var_data[3].bits[1].kind == PC.FIXED1
    @test postsolver.var_data[3].bits[1].ref_var_id === nothing
    @test postsolver.var_data[1].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[1].bits[1].ref_var_id == 3
    @test postsolver.var_data[2].bits[1].kind == PC.NEGATED_BINVAR
    @test postsolver.var_data[2].bits[1].ref_var_id == 3
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 2.0, 2 => 3.0)) ==
          Dict{PC.VarId, Float64}(1 => 5.0, 2 => 6.0)
end

@testset "postsolve recursively resolves nested helper vars introduced by repeated parity patterns" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(0.0, 5.0),
        4 => PC.IntVar(0.0, 3.0),
    )
    model = PC.QPModel(vars, PC.Constraint[], parity_empty_objective(), :min)
    manager = PC.PropagationManager([1, 2, 4])
    postsolver = PC.ParityPostsolver(keys(vars))

    PC.add_equivalence!(manager, PC.VarLit(1, false), PC.VarLit(2, false))
    PC.update_sccs!(manager)
    while PC.pop_substitution!(manager) !== nothing
    end

    @test PC.fix_parity_patterns!(model, manager, postsolver) == 2
    @test haskey(model.vars, 5)
    @test haskey(postsolver.var_data, 5)
    @test postsolver.var_data[1].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[1].bits[1].ref_var_id == 5
    @test postsolver.var_data[2].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[2].bits[1].ref_var_id == 5

    PC.add_equivalence!(manager, PC.VarLit(5, false), PC.VarLit(4, false))
    PC.update_sccs!(manager)
    while PC.pop_substitution!(manager) !== nothing
    end

    @test PC.fix_parity_patterns!(model, manager, postsolver) == 2
    @test haskey(model.vars, 6)
    PC.fix_vars!(model, postsolver)

    @test !haskey(model.vars, 5)
    @test postsolver.var_data[5].current_var_id === nothing
    @test postsolver.var_data[5].fixed_high_order == 0.0
    @test length(postsolver.var_data[5].bits) == 1
    @test postsolver.var_data[5].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[5].bits[1].ref_var_id == 6
    @test postsolver.var_data[1].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[1].bits[1].ref_var_id == 5
    @test postsolver.var_data[2].bits[1].kind == PC.BINVAR
    @test postsolver.var_data[2].bits[1].ref_var_id == 5
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 1.0, 2 => 2.0, 4 => 1.0, 6 => 1.0)) ==
          Dict{PC.VarId, Float64}(1 => 3.0, 2 => 5.0, 4 => 3.0)
end

@testset "parity_presolve! records repeated parity reductions in bit order" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 4.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(1.0, 1, 1)], ParityLinTerm[]),
        4.0,
        4.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    postsolver = PC.ParityPostsolver(keys(vars))

    PC.parity_presolve!(model, postsolver)

    var_data = postsolver.var_data[1]
    @test !model.infeasible
    @test isempty(model.vars)
    @test var_data.current_var_id === nothing
    @test var_data.fixed_high_order == 0.0
    @test length(var_data.bits) == 2
    @test var_data.bits[1].kind == PC.FIXED0
    @test var_data.bits[2].kind == PC.FIXED1
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}()) == Dict{PC.VarId, Float64}(1 => 2.0)
end

@testset "parity_presolve! reconstructs shifted binaries normalized before parity" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(3.0, 4.0))
    model = PC.QPModel(vars, PC.Constraint[], parity_empty_objective(), :min)
    postsolver = PC.ParityPostsolver(keys(vars))

    stats = PC.parity_presolve!(model, postsolver)

    @test !model.infeasible
    @test stats.changed
    @test stats.domains_changed
    @test isempty(stats.coefficient_changed_constraint_ids)
    @test stats.fixed_parities == 0
    @test stats.pattern_rewritten_vars == 0
    @test !stats.infeasible
    @test model.vars[1] == PC.IntVar(0.0, 1.0)
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 0.0)) == Dict{PC.VarId, Float64}(1 => 3.0)
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}(1 => 1.0)) == Dict{PC.VarId, Float64}(1 => 4.0)
end

@testset "parity_presolve! reconstructs shifted binaries removed during normalization" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(3.0, 4.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1)]),
        4.0,
        4.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    postsolver = PC.ParityPostsolver(keys(vars))

    PC.parity_presolve!(model, postsolver)

    @test !model.infeasible
    @test isempty(model.vars)
    @test isempty(model.cons)
    @test PC.postsolve(postsolver, Dict{PC.VarId, Float64}()) == Dict{PC.VarId, Float64}(1 => 4.0)
end

@testset "parity_presolve_phase! normalizes gcd-reducible constraints at phase entry" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0))
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(2.0, 1)]),
        2.0,
        2.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)
    propagator = PC.PropagationManager(PC.VarId[])

    stats = PC.parity_presolve_phase!(model, propagator)

    @test !model.infeasible
    @test isempty(model.vars)
    @test isempty(model.cons)
    @test stats.changed
end

@testset "parity_presolve! accepts shifted binary inputs without pre-normalization" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(2.0, 3.0),
        2 => PC.IntVar(0.0, 1.0),
    )
    con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2)]),
        3.0,
        3.0,
    )
    model = PC.QPModel(vars, [con], parity_empty_objective(), :min)

    PC.parity_presolve!(model)

    @test !model.infeasible
    @test all(PC.is_binary(var) for var in values(model.vars))
end

@testset "parity_presolve! supports parity strategy presets" begin
    mod2_model = parity_single_distance_model()
    mod2_result = PC.parity_presolve!(
        mod2_model;
        collect_stats = true,
        parity_strategy = :mod2_basic,
    )
    @test !mod2_model.infeasible
    @test !mod2_result.infeasible
    @test mod2_result.parity_stats.num_xorand_constraints_generated == 0

    mod4_basic_model = parity_single_distance_model()
    mod4_basic_result = PC.parity_presolve!(
        mod4_basic_model;
        collect_stats = true,
        parity_strategy = :mod4_basic,
    )
    @test !mod4_basic_model.infeasible
    @test !mod4_basic_result.infeasible
    @test mod4_basic_result.parity_stats.num_xorand_constraints_generated > 0

    full_model = parity_single_distance_model()
    full_result = PC.parity_presolve!(
        full_model;
        collect_stats = true,
        parity_strategy = :full,
    )
    @test !full_model.infeasible
    @test !full_result.infeasible
    @test full_result.parity_stats.num_xorand_constraints_generated > 0

    @test_throws ArgumentError PC.parity_presolve!(
        parity_single_distance_model();
        parity_strategy = :unknown,
    )
end
