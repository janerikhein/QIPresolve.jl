using Test
using Graphs: has_edge, ne
import QIPresolve.PresolvingCore as PC

function parity_identity_maps(n::Int)
    pos_to_var = collect(1:n)
    var_to_pos = Dict(i => i for i in 1:n)
    return pos_to_var, var_to_pos
end

function parity_symmetric_bitmatrix(n::Int, edges::Vector{Tuple{Int, Int}})
    mat = falses(n, n)
    for (i, j) in edges
        mat[i, j] = true
        mat[j, i] = true
    end
    return mat
end

function parity_scc_pos(manager::PC.PropagationManager, lit::PC.VarLit)
    rep = PC.repr!(manager, manager.lit_to_pos[lit])
    return manager.rep_pos_to_scc_pos[rep]
end

@testset "propagate! drains stale fixing events for absent parity variables" begin
    model = PC.ParityModel(
        Dict{PC.VarId, Int}(1 => 1),
        PC.VarId[1],
        [PC.XorConstraint(BitVector([1]), false)],
    )
    manager = PC.PropagationManager(PC.VarId[1, 2])
    PC.fix_var!(manager, 2, true)

    @test !isempty(manager.pending_fixings)
    @test PC.propagate!(model, manager) === model
    @test isempty(manager.pending_fixings)
    @test !model.infeasible
end

@testset "propagate! drains stale substitution events for absent parity variables" begin
    absent_source_model = PC.ParityModel(
        Dict{PC.VarId, Int}(1 => 1),
        PC.VarId[1],
        [PC.XorConstraint(BitVector([1]), false)],
    )
    absent_source_manager = PC.PropagationManager(PC.VarId[1, 2, 3])
    PC.add_equivalence!(
        absent_source_manager,
        PC.VarLit(2, false),
        PC.VarLit(3, false),
    )

    @test PC.propagate!(absent_source_model, absent_source_manager) === absent_source_model
    @test isempty(absent_source_manager.pending_substitutions)
    @test !absent_source_model.infeasible

    absent_target_model = PC.ParityModel(
        Dict{PC.VarId, Int}(2 => 1),
        PC.VarId[2],
        [PC.XorConstraint(BitVector([1]), false)],
    )
    absent_target_manager = PC.PropagationManager(PC.VarId[1, 2])
    PC.add_equivalence!(
        absent_target_manager,
        PC.VarLit(1, false),
        PC.VarLit(2, false),
    )

    @test PC.propagate!(absent_target_model, absent_target_manager) === absent_target_model
    @test isempty(absent_target_manager.pending_substitutions)
    @test !absent_target_model.infeasible
end

@testset "gauss_jordan_xor! only pivots and eliminates pure XOR rows" begin
    pos_to_var, var_to_pos = parity_identity_maps(3)
    xor1 = PC.XorConstraint(BitVector([1, 1, 0]), false)
    xor2 = PC.XorConstraint(BitVector([1, 0, 1]), true)
    xor_and = PC.XorConstraint(
        BitVector([1, 0, 0]),
        parity_symmetric_bitmatrix(3, [(2, 3)]),
        true,
    )

    model = PC.ParityModel(var_to_pos, pos_to_var, [xor1, xor2, deepcopy(xor_and)])
    PC.gauss_jordan_xor!(model)

    @test model.pivots[1] == (1, nothing)
    @test model.pivots[2] == (2, nothing)
    @test model.pivots[3] === nothing

    @test model.cons[1].par == BitVector([1, 0, 1])
    @test model.cons[1].rhs
    @test model.cons[2].par == BitVector([0, 1, 1])
    @test model.cons[2].rhs
    @test model.cons[3].par == xor_and.par
    @test model.cons[3].conj == xor_and.conj
    @test model.cons[3].rhs == xor_and.rhs
end

@testset "gauss_jordan_xor_and! only pivots and eliminates XOR-AND rows" begin
    pos_to_var, var_to_pos = parity_identity_maps(4)
    xor_con = PC.XorConstraint(BitVector([1, 1, 0, 0]), false)
    xor_and1 = PC.XorConstraint(
        BitVector([1, 0, 0, 0]),
        parity_symmetric_bitmatrix(4, [(1, 2)]),
        false,
    )
    xor_and2 = PC.XorConstraint(
        BitVector([0, 1, 0, 0]),
        parity_symmetric_bitmatrix(4, [(1, 2), (2, 3)]),
        true,
    )

    model = PC.ParityModel(var_to_pos, pos_to_var, [deepcopy(xor_con), xor_and1, xor_and2])
    PC.gauss_jordan_xor_and!(model)

    @test model.pivots[1] === nothing
    @test model.pivots[2] == (2, 1)
    @test model.pivots[3] == (3, 2)

    @test model.cons[1].par == xor_con.par
    @test model.cons[1].rhs == xor_con.rhs
    @test model.cons[2].conj == parity_symmetric_bitmatrix(4, [(1, 2)])
    @test model.cons[3].conj == parity_symmetric_bitmatrix(4, [(2, 3)])
    @test model.cons[3].par == BitVector([1, 1, 0, 0])
    @test model.cons[3].rhs
end

@testset "substitute_parity_pivots! clears stored parity pivots from non-owner rows" begin
    pos_to_var, var_to_pos = parity_identity_maps(4)
    xor_pivot = PC.XorConstraint(BitVector([1, 1, 0, 0]), false)
    xor_and_pivot = PC.XorConstraint(
        BitVector([1, 0, 0, 0]),
        parity_symmetric_bitmatrix(4, [(2, 3)]),
        false,
    )
    xor_unpivoted = PC.XorConstraint(BitVector([1, 0, 1, 0]), true)
    xor_and_unpivoted = PC.XorConstraint(
        BitVector([1, 0, 0, 0]),
        parity_symmetric_bitmatrix(4, [(3, 4)]),
        true,
    )

    model = PC.ParityModel(
        var_to_pos,
        pos_to_var,
        [xor_pivot, xor_and_pivot, xor_unpivoted, xor_and_unpivoted],
    )
    model.pivots[1] = (1, nothing)
    model.pivots[2] = (2, 3)

    PC.substitute_parity_pivots!(model)

    @test model.cons[2].par == BitVector([0, 1, 0, 0])
    @test model.cons[2].conj == parity_symmetric_bitmatrix(4, [(2, 3)])
    @test !model.cons[2].rhs

    @test model.cons[3].par == BitVector([0, 1, 1, 0])
    @test model.cons[3].rhs

    @test model.cons[4].par == BitVector([0, 1, 0, 0])
    @test model.cons[4].conj == parity_symmetric_bitmatrix(4, [(3, 4)])
    @test model.cons[4].rhs
end

@testset "substitute_parity_pivots! reduces same-type pivot rows retaining a parity pivot" begin
    pos_to_var, var_to_pos = parity_identity_maps(3)
    xor_pivot1 = PC.XorConstraint(BitVector([1, 1, 0]), false)
    xor_pivot2 = PC.XorConstraint(BitVector([1, 0, 1]), true)

    model = PC.ParityModel(var_to_pos, pos_to_var, [xor_pivot1, xor_pivot2])
    model.pivots[1] = (1, nothing)
    model.pivots[2] = (3, nothing)

    PC.substitute_parity_pivots!(model)

    @test model.cons[2].par == BitVector([0, 1, 1])
    @test model.cons[2].rhs
    @test model.pivots[2] == (3, nothing)
    @test !model.infeasible
end

@testset "substitute_parity_pivots! ignores conjunctive pivot owners" begin
    pos_to_var, var_to_pos = parity_identity_maps(3)
    xor_and_pivot = PC.XorConstraint(
        BitVector([1, 0, 0]),
        parity_symmetric_bitmatrix(3, [(1, 2)]),
        false,
    )
    xor_row = PC.XorConstraint(BitVector([1, 1, 0]), true)

    model = PC.ParityModel(var_to_pos, pos_to_var, [deepcopy(xor_and_pivot), deepcopy(xor_row)])
    model.pivots[1] = (1, 2)

    PC.substitute_parity_pivots!(model)

    @test model.cons[1].par == xor_and_pivot.par
    @test model.cons[1].conj == xor_and_pivot.conj
    @test model.cons[1].rhs == xor_and_pivot.rhs
    @test model.cons[2].par == xor_row.par
    @test model.cons[2].rhs == xor_row.rhs
end

@testset "substitute_pivots_in_conjunctive_terms! rewrites non-owner conjunction terms only" begin
    pos_to_var, var_to_pos = parity_identity_maps(4)
    xor_pivot = PC.XorConstraint(BitVector([1, 1, 1, 0]), true)
    target = PC.XorConstraint(
        falses(4),
        parity_symmetric_bitmatrix(4, [(1, 4)]),
        false,
    )
    linear_only = PC.XorConstraint(BitVector([1, 0, 0, 1]), false)

    model = PC.ParityModel(var_to_pos, pos_to_var, [deepcopy(xor_pivot), target, deepcopy(linear_only)])
    model.pivots[1] = (1, nothing)

    out = PC.substitute_pivots_in_conjunctive_terms!(model)

    @test out === model
    @test model.cons[1].par == xor_pivot.par
    @test model.cons[1].rhs == xor_pivot.rhs
    @test model.cons[2].par == BitVector([0, 0, 0, 1])
    @test model.cons[2].conj == parity_symmetric_bitmatrix(4, [(2, 4), (3, 4)])
    @test !model.cons[2].rhs
    @test model.cons[3].par == linear_only.par
    @test model.cons[3].rhs == linear_only.rhs
end

@testset "substitute_pivots_in_conjunctive_terms! revalidates changed pivots and composes multiple substitutions" begin
    pos_to_var, var_to_pos = parity_identity_maps(5)
    xor_pivot1 = PC.XorConstraint(BitVector([1, 1, 0, 0, 0]), true)
    xor_pivot2 = PC.XorConstraint(BitVector([0, 0, 1, 0, 1]), false)
    target = PC.XorConstraint(
        falses(5),
        parity_symmetric_bitmatrix(5, [(1, 4), (3, 4)]),
        false,
    )

    model = PC.ParityModel(var_to_pos, pos_to_var, [xor_pivot1, xor_pivot2, target])
    model.pivots[1] = (1, nothing)
    model.pivots[2] = (3, nothing)
    model.pivots[3] = (1, 4)

    PC.substitute_pivots_in_conjunctive_terms!(model)

    @test model.cons[3].par == BitVector([0, 0, 0, 1, 0])
    @test model.cons[3].conj == parity_symmetric_bitmatrix(5, [(2, 4), (4, 5)])
    @test !model.cons[3].rhs
    @test model.pivots[3] === nothing
end

@testset "propagate! rewrites full bipartite rows into two XOR constraints" begin
    pos_to_var, var_to_pos = parity_identity_maps(7)
    bipartite = PC.XorConstraint(
        falses(7),
        parity_symmetric_bitmatrix(7, [(1, 4), (1, 5), (1, 6), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6)]),
        true,
    )
    unsplittable = PC.XorConstraint(
        BitVector([1, 0, 0, 0, 0, 1, 1]),
        false,
    )
    untouched = PC.XorConstraint(
        BitVector([0, 1, 0, 1, 0, 0, 1]),
        false,
    )

    model = PC.ParityModel(var_to_pos, pos_to_var, [unsplittable, bipartite, untouched])
    model.pivots[1] = (1, nothing)
    model.pivots[2] = (4, 1)
    model.pivots[3] = (4, nothing)
    manager = PC.PropagationManager(model.pos_to_var_id)

    out = PC.propagate!(model, manager)

    @test out === model
    @test length(model.cons) == 4
    @test model.pivots == PC.PivotSlot[(1, nothing), nothing, nothing, (4, nothing)]

    @test model.cons[1].par == unsplittable.par
    @test model.cons[1].rhs == unsplittable.rhs

    @test model.cons[2].meta.is_pure_xor
    @test model.cons[3].meta.is_pure_xor
    @test model.cons[2].rhs
    @test model.cons[3].rhs

    expected_split = Set([BitVector([1, 1, 1, 0, 0, 0, 0]), BitVector([0, 0, 0, 1, 1, 1, 0])])
    @test Set([model.cons[2].par, model.cons[3].par]) == expected_split

    @test model.cons[4].par == untouched.par
    @test model.cons[4].rhs == untouched.rhs
    @test all(!con.meta.requires_prop for con in model.cons)
end

@testset "propagate! rewrites tripartite rows with nonempty Z into XOR constraints" begin
    pos_to_var, var_to_pos = parity_identity_maps(4)
    tripartite = PC.XorConstraint(
        BitVector([1, 0, 0, 0]),
        parity_symmetric_bitmatrix(4, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)]),
        true,
    )
    model = PC.ParityModel(var_to_pos, pos_to_var, [tripartite])
    manager = PC.PropagationManager(model.pos_to_var_id)

    PC.propagate!(model, manager)

    @test length(model.cons) == 2
    @test model.pivots == PC.PivotSlot[nothing, nothing]
    @test model.cons[1].meta.is_pure_xor
    @test model.cons[1].par == BitVector([1, 0, 1, 1])
    @test model.cons[1].rhs
    @test model.cons[2].meta.is_pure_xor
    @test model.cons[2].par == BitVector([0, 1, 1, 1])
    @test !model.cons[2].rhs
    @test all(!con.meta.requires_prop for con in model.cons)

    basic_tripartite = PC.XorConstraint(
        BitVector([1, 0, 0, 0]),
        parity_symmetric_bitmatrix(4, [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)]),
        true,
    )
    basic_model = PC.ParityModel(var_to_pos, pos_to_var, [basic_tripartite])
    basic_manager = PC.PropagationManager(basic_model.pos_to_var_id)

    PC.propagate!(basic_model, basic_manager; parity_strategy = :mod4_basic)

    @test length(basic_model.cons) == 1
    @test basic_model.pivots == PC.PivotSlot[nothing]
    @test basic_model.cons[1].par == basic_tripartite.par
    @test basic_model.cons[1].conj == basic_tripartite.conj
    @test basic_model.cons[1].rhs == basic_tripartite.rhs
    @test !basic_model.cons[1].meta.requires_prop
    @test ne(basic_manager.graph) == 0
end

@testset "propagate! adds implication-only tripartite matches without removing the row" begin
    pos_to_var, var_to_pos = parity_identity_maps(2)
    implication_row = PC.XorConstraint(
        BitVector([0, 1]),
        parity_symmetric_bitmatrix(2, [(1, 2)]),
        false,
    )
    model = PC.ParityModel(var_to_pos, pos_to_var, [implication_row])
    manager = PC.PropagationManager(model.pos_to_var_id)

    PC.propagate!(model, manager)

    @test length(model.cons) == 1
    @test model.cons[1].par == implication_row.par
    @test model.cons[1].conj == implication_row.conj
    @test model.cons[1].rhs == implication_row.rhs
    @test !model.cons[1].meta.requires_prop
    @test !model.infeasible

    x_false = parity_scc_pos(manager, PC.VarLit(1, true))
    y_false = parity_scc_pos(manager, PC.VarLit(2, true))
    y_true = parity_scc_pos(manager, PC.VarLit(2, false))
    x_true = parity_scc_pos(manager, PC.VarLit(1, false))

    @test has_edge(manager.graph, x_false, y_false)
    @test has_edge(manager.graph, y_true, x_true)
    @test PC.fixed_value(manager, 1) === nothing
    @test PC.fixed_value(manager, 2) === nothing
end

@testset "propagate! rewrites one tripartite side into a fixing and the other into an XOR row" begin
    pos_to_var, var_to_pos = parity_identity_maps(4)
    mixed_row = PC.XorConstraint(
        BitVector([0, 1, 1, 1]),
        parity_symmetric_bitmatrix(4, [(1, 2), (1, 3), (1, 4)]),
        true,
    )
    model = PC.ParityModel(var_to_pos, pos_to_var, [mixed_row])
    manager = PC.PropagationManager(model.pos_to_var_id)

    PC.propagate!(model, manager)

    @test length(model.cons) == 1
    @test model.pivots == PC.PivotSlot[nothing]
    @test model.cons[1].meta.is_pure_xor
    @test model.cons[1].par == BitVector([0, 1, 1, 1])
    @test model.cons[1].rhs
    @test PC.fixed_value(manager, 1) == false
    @test !model.infeasible
end

@testset "propagate! rewrites both tripartite sides into direct fixings" begin
    pos_to_var, var_to_pos = parity_identity_maps(2)
    fixing_row = PC.XorConstraint(
        BitVector([1, 1]),
        parity_symmetric_bitmatrix(2, [(1, 2)]),
        false,
    )
    model = PC.ParityModel(var_to_pos, pos_to_var, [fixing_row])
    manager = PC.PropagationManager(model.pos_to_var_id)

    PC.propagate!(model, manager)

    @test isempty(model.cons)
    @test isempty(model.pivots)
    @test PC.fixed_value(manager, 1) == false
    @test PC.fixed_value(manager, 2) == false
    @test !model.infeasible
end

@testset "propagate! leaves bad tripartite linear supports unchanged" begin
    pos_to_var, var_to_pos = parity_identity_maps(3)
    unchanged = PC.XorConstraint(
        BitVector([1, 1, 0]),
        parity_symmetric_bitmatrix(3, [(1, 2), (1, 3), (2, 3)]),
        true,
    )
    model = PC.ParityModel(var_to_pos, pos_to_var, [deepcopy(unchanged)])
    model.pivots[1] = (1, 2)
    manager = PC.PropagationManager(model.pos_to_var_id)

    PC.propagate!(model, manager)

    @test length(model.cons) == 1
    @test model.pivots == PC.PivotSlot[(1, 2)]
    @test model.cons[1].par == unchanged.par
    @test model.cons[1].conj == unchanged.conj
    @test model.cons[1].rhs == unchanged.rhs
    @test !model.cons[1].meta.requires_prop
    @test !model.infeasible
end

@testset "cleanup! removes empty rows, tracks infeasibility, and keeps pivots aligned" begin
    pos_to_var, var_to_pos = parity_identity_maps(4)
    empty_feasible = PC.XorConstraint(falses(4), false)
    keep1 = PC.XorConstraint(BitVector([1, 0, 0, 1]), true)
    empty_infeasible = PC.XorConstraint(falses(4), true)
    keep2 = PC.XorConstraint(
        BitVector([0, 1, 0, 0]),
        parity_symmetric_bitmatrix(4, [(2, 3)]),
        false,
    )
    stale_empty = PC.XorConstraint(
        BitVector([1, 0, 0, 0]),
        parity_symmetric_bitmatrix(4, [(1, 2)]),
        true,
    )
    PC.xor_con!(stale_empty, stale_empty)

    model = PC.ParityModel(
        var_to_pos,
        pos_to_var,
        [empty_feasible, keep1, empty_infeasible, keep2, stale_empty],
    )
    model.pivots[1] = (1, nothing)
    model.pivots[2] = (4, nothing)
    model.pivots[3] = (2, nothing)
    model.pivots[4] = (3, 2)
    model.pivots[5] = (1, 2)

    out = PC.cleanup!(model)

    @test out === model
    @test model.infeasible
    @test length(model.cons) == 2
    @test model.pivots == PC.PivotSlot[(4, nothing), (3, 2)]

    @test model.cons[1].par == keep1.par
    @test model.cons[1].rhs == keep1.rhs
    @test model.cons[2].par == keep2.par
    @test model.cons[2].conj == keep2.conj
    @test model.cons[2].rhs == keep2.rhs
end

@testset "fix_var!(model) applies to pure XOR and XOR-AND rows" begin
    pos_to_var, var_to_pos = parity_identity_maps(3)
    xor_con = PC.XorConstraint(BitVector([1, 0, 1]), true)
    xor_and = PC.XorConstraint(
        BitVector([0, 1, 0]),
        parity_symmetric_bitmatrix(3, [(1, 2)]),
        false,
    )

    model = PC.ParityModel(var_to_pos, pos_to_var, [xor_con, xor_and])
    out = PC.fix_var!(model, 1, true)

    @test out === model
    @test model.cons[1].par == BitVector([0, 0, 1])
    @test !model.cons[1].rhs
    @test model.cons[1].meta.requires_update
    @test model.cons[1].meta.requires_prop

    @test model.cons[2].par == falses(3)
    @test model.cons[2].conj == parity_symmetric_bitmatrix(3, Tuple{Int, Int}[])
    @test !model.cons[2].rhs
    @test model.cons[2].meta.requires_update
    @test model.cons[2].meta.requires_prop

    PC.ensure_updated!(model.cons[2])
    @test model.cons[2].conj === nothing
end

@testset "propagate! applies substitutions from two-term XOR rows" begin
    pos_to_var, var_to_pos = parity_identity_maps(4)
    eq_row = PC.XorConstraint(BitVector([1, 1, 0, 0]), true)
    other_row = PC.XorConstraint(BitVector([0, 1, 1, 1]), false)

    model = PC.ParityModel(var_to_pos, pos_to_var, [eq_row, other_row])
    manager = PC.PropagationManager(model.pos_to_var_id)

    out = PC.propagate!(model, manager)

    @test out === model
    @test length(model.cons) == 1
    PC.ensure_updated!(model.cons[1])
    @test !model.infeasible
    @test model.cons[1].par == BitVector([1, 0, 1, 1])
    @test model.cons[1].rhs
    @test all(!con.meta.requires_prop for con in model.cons)
end

@testset "propagate! applies implications from mixed two-term rows" begin
    pos_to_var, var_to_pos = parity_identity_maps(3)
    mixed_row = PC.XorConstraint(
        BitVector([0, 0, 1]),
        parity_symmetric_bitmatrix(3, [(1, 2)]),
        false,
    )
    fix_row = PC.XorConstraint(BitVector([0, 0, 1]), true)

    model = PC.ParityModel(var_to_pos, pos_to_var, [mixed_row, fix_row])
    manager = PC.PropagationManager(model.pos_to_var_id)

    out = PC.propagate!(model, manager)

    @test out === model
    @test !model.infeasible
    @test PC.fixed_value(manager, 1) == true
    @test PC.fixed_value(manager, 2) == true
    @test PC.fixed_value(manager, 3) == true
    @test isempty(model.cons)
    @test isempty(model.pivots)
end

@testset "propagate! re-enqueues fixings across passes when substitutions create fresh occurrences" begin
    pos_to_var, var_to_pos = parity_identity_maps(3)
    fix_row = PC.XorConstraint(BitVector([1, 0, 0]), false)
    eq_row = PC.XorConstraint(BitVector([1, 1, 0]), false)
    mixed_row = PC.XorConstraint(
        BitVector([0, 1, 0]),
        parity_symmetric_bitmatrix(3, [(2, 3)]),
        false,
    )

    model = PC.ParityModel(var_to_pos, pos_to_var, [fix_row, eq_row, mixed_row])
    manager = PC.PropagationManager(model.pos_to_var_id)

    out = PC.propagate!(model, manager)

    @test out === model
    @test !model.infeasible
    @test isempty(model.cons)
    @test isempty(model.pivots)
end

@testset "propagate! does not make the p11/p1 fragment infeasible" begin
    pos_to_var = [1, 7, 11, 21]
    var_to_pos = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var))
    conj = falses(4, 4)
    conj[1, 4] = true
    conj[4, 1] = true

    model = PC.ParityModel(
        var_to_pos,
        pos_to_var,
        [
            PC.XorConstraint(BitVector([1, 0, 0, 1]), false),
            PC.XorConstraint(falses(4), conj, false),
            PC.XorConstraint(BitVector([1, 1, 0, 0]), true),
            PC.XorConstraint(BitVector([0, 1, 0, 0]), true),
            PC.XorConstraint(BitVector([1, 0, 1, 0]), true),
        ],
    )
    manager = PC.PropagationManager(model.pos_to_var_id)

    PC.propagate!(model, manager)

    @test !model.infeasible
    @test all(begin
        PC.ensure_updated!(con)
        !(PC.constraint_nnz(con) == 0 && con.rhs)
    end for con in model.cons)
    @test all(!con.meta.requires_prop for con in model.cons)
end

@testset "propagate! reaches a fixpoint and revalidates changed pivots" begin
    pos_to_var, var_to_pos = parity_identity_maps(5)
    eq_row = PC.XorConstraint(BitVector([1, 1, 0, 0, 0]), false)
    fix_row = PC.XorConstraint(BitVector([0, 1, 0, 0, 0]), true)
    untouched = PC.XorConstraint(BitVector([0, 0, 1, 1, 1]), true)

    model = PC.ParityModel(var_to_pos, pos_to_var, [eq_row, fix_row, untouched])
    model.pivots[1] = (1, nothing)
    model.pivots[3] = (3, nothing)
    manager = PC.PropagationManager(model.pos_to_var_id)

    PC.propagate!(model, manager)

    @test length(model.cons) == 1
    PC.ensure_updated!(model.cons[1])

    @test !model.infeasible
    @test model.cons[1].par == untouched.par
    @test model.cons[1].rhs == untouched.rhs
    @test model.pivots == PC.PivotSlot[(3, nothing)]
    @test all(!con.meta.requires_prop for con in model.cons)
end

@testset "show(::ParityModel) prints grouped algebraic constraints with model-owned variable ids" begin
    pos_to_var = [20, 40, 30, 10]
    var_to_pos = Dict(var_id => pos for (pos, var_id) in enumerate(pos_to_var))

    xor1 = PC.XorConstraint(BitVector([1, 1, 0, 0]), false)
    xor_and = PC.XorConstraint(
        BitVector([0, 1, 0, 0]),
        parity_symmetric_bitmatrix(4, [(1, 3)]),
        true,
    )
    empty_xor = PC.XorConstraint(falses(4), false)
    xor2 = PC.XorConstraint(BitVector([0, 0, 0, 1]), true)

    model = PC.ParityModel(var_to_pos, pos_to_var, [xor_and, xor1, empty_xor, xor2])
    model.infeasible = true
    model.pivots[1] = (1, 3)
    model.pivots[4] = (4, nothing)

    shown = repr("text/plain", model)

    @test occursin("ParityModel", shown)
    @test occursin("Variables: 4", shown)
    @test occursin("Constraints: 4", shown)
    @test occursin("XOR constraints: 3", shown)
    @test occursin("XOR-AND constraints: 1", shown)
    @test occursin("Infeasible: true", shown)
    @test occursin("Pivoted rows: 2", shown)
    @test occursin("XOR pivots: 1", shown)
    @test occursin("XOR-AND pivots: 1", shown)
    @test occursin("p20 ⊕ p40 = 0", shown)
    @test occursin("p10 = 1", shown)
    @test occursin("0 = 0", shown)
    @test occursin("p40 ⊕ (p20 ∧ p30) = 1", shown)

    xor_section_idx = findfirst("XOR constraints (3):", shown)
    xor1_idx = findfirst("p20 ⊕ p40 = 0", shown)
    empty_idx = findfirst("0 = 0", shown)
    xor2_idx = findfirst("p10 = 1", shown)
    xor_and_section_idx = findfirst("XOR-AND constraints (1):", shown)
    xor_and_idx = findfirst("p40 ⊕ (p20 ∧ p30) = 1", shown)

    @test xor_section_idx !== nothing
    @test xor1_idx !== nothing
    @test empty_idx !== nothing
    @test xor2_idx !== nothing
    @test xor_and_section_idx !== nothing
    @test xor_and_idx !== nothing
    @test xor_section_idx < xor1_idx < empty_idx < xor2_idx < xor_and_section_idx < xor_and_idx
end
