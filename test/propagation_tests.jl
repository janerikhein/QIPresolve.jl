using Test
using Graphs: has_edge, inneighbors, ne, nv, outneighbors, topological_sort
import QIPresolve.PresolvingCore as PC

function prepare_substitution_manager()
    manager = PC.PropagationManager([1, 2, 3])
    PC.add_equivalence!(manager, PC.VarLit(1, false), PC.VarLit(2, false))
    PC.add_implication!(manager, PC.VarLit(1, false), PC.VarLit(3, true))
    PC.add_implication!(manager, PC.VarLit(3, false), PC.VarLit(1, false))
    PC.update_sccs!(manager)
    drain_substitutions!(manager)
    return manager
end

function drain_fixings!(manager)
    items = Tuple{Int, Bool}[]
    while true
        item = PC.pop_fixing!(manager)
        item === nothing && return items
        push!(items, item)
    end
end

function drain_substitutions!(manager)
    items = Tuple{Int, Int, Bool}[]
    while true
        item = PC.pop_substitution!(manager)
        item === nothing && return items
        push!(items, item)
    end
end

function propagation_edge_matrix(n::Int, edges::Vector{Tuple{Int, Int}})
    mat = falses(n, n)
    for (i, j) in edges
        mat[i, j] = true
        mat[j, i] = true
    end
    return mat
end

@testset "VarLit helpers and constructor" begin
    lit = PC.VarLit(7, false)
    @test PC.negated(lit) == PC.VarLit(7, true)
    @test PC.negated(PC.negated(lit)) == lit

    manager = PC.PropagationManager([10, 20, 30])

    @test manager.nvars == 3
    @test manager.nreps == 6
    @test nv(manager.graph) == 6
    @test ne(manager.graph) == 0
    @test manager.parent_pos == collect(1:6)
    @test manager.scc_pos_to_rep_pos == collect(1:6)
    @test manager.rep_pos_to_scc_pos == Dict(i => i for i in 1:6)
    @test manager.lit_labels == fill(PC.UNDEF, 6)
    @test PC.pop_fixing!(manager) === nothing
    @test PC.pop_substitution!(manager) === nothing

    @test manager.lit_to_pos[PC.VarLit(10, false)] == 1
    @test manager.lit_to_pos[PC.VarLit(10, true)] == 2
    @test manager.lit_to_pos[PC.VarLit(20, false)] == 3
    @test manager.lit_to_pos[PC.VarLit(20, true)] == 4
    @test manager.lit_to_pos[PC.VarLit(30, false)] == 5
    @test manager.lit_to_pos[PC.VarLit(30, true)] == 6
    @test manager.pos_to_lit == [
        PC.VarLit(10, false),
        PC.VarLit(10, true),
        PC.VarLit(20, false),
        PC.VarLit(20, true),
        PC.VarLit(30, false),
        PC.VarLit(30, true),
    ]

    @test_throws AssertionError PC.PropagationManager([1, 1])
end

@testset "repr! and queue helpers" begin
    manager = PC.PropagationManager([1, 2, 3])
    manager.parent_pos[6] = 4
    manager.parent_pos[5] = 2
    manager.parent_pos[4] = 2
    manager.parent_pos[2] = 2

    @test PC.repr!(manager, 6) == 2
    @test manager.parent_pos[6] == 2
    @test manager.parent_pos[4] == 2

    @test PC.enqueue_fixing!(manager, 1, true)
    @test !PC.enqueue_fixing!(manager, 1, true)
    @test PC.enqueue_fixing!(manager, 2, false)
    @test PC.pop_fixing!(manager) == (1, true)
    @test PC.pop_fixing!(manager) == (2, false)
    @test PC.pop_fixing!(manager) === nothing

    @test PC.enqueue_substitution!(manager, 1)
    @test !PC.enqueue_substitution!(manager, 1)
    @test PC.enqueue_substitution!(manager, 3)
    @test PC.pop_substitution!(manager) == (1, 1, false)
    @test PC.pop_substitution!(manager) == (3, 1, true)
    @test PC.pop_substitution!(manager) === nothing
end

@testset "ensure_literals! appends fresh singleton literals only for missing vars" begin
    manager = PC.PropagationManager([1])

    @test PC.ensure_literals!(manager, [1, 2, 3]) === manager
    @test manager.nvars == 3
    @test manager.nreps == 6
    @test nv(manager.graph) == 6
    @test ne(manager.graph) == 0
    @test manager.lit_to_pos[PC.VarLit(1, false)] == 1
    @test manager.lit_to_pos[PC.VarLit(1, true)] == 2
    @test manager.lit_to_pos[PC.VarLit(2, false)] == 3
    @test manager.lit_to_pos[PC.VarLit(2, true)] == 4
    @test manager.lit_to_pos[PC.VarLit(3, false)] == 5
    @test manager.lit_to_pos[PC.VarLit(3, true)] == 6

    for pos in 1:6
        @test manager.parent_pos[pos] == pos
        @test manager.rep_pos_to_scc_pos[pos] == pos
        @test manager.scc_pos_to_rep_pos[pos] == pos
        @test manager.lit_labels[pos] == PC.UNDEF
    end
end

@testset "Label helpers" begin
    manager = PC.PropagationManager([1])

    @test PC.set_lit_label!(manager, 1, PC.TRUE)
    @test manager.lit_labels[1] == PC.TRUE
    @test PC.pop_fixing!(manager) === nothing

    @test !PC.set_lit_label!(manager, 1, PC.TRUE)
    @test PC.set_lit_label!(manager, 2, PC.FALSE)
    @test manager.lit_labels[2] == PC.FALSE
    @test PC.pop_fixing!(manager) == (1, true)
    @test PC.pop_fixing!(manager) === nothing

    @test PC.set_lit_label!(manager, 1, PC.FALSE)
    @test manager.lit_labels[1] == PC.FALSE

    manager2 = PC.PropagationManager([2])
    manager2.lit_labels[1] = PC.FALSE
    manager2.lit_labels[2] = PC.TRUE
    PC.maybe_enqueue_fixing!(manager2, 1)
    @test PC.pop_fixing!(manager2) == (2, false)
    PC.maybe_enqueue_fixing!(manager2, 2)
    @test PC.pop_fixing!(manager2) === nothing
end

@testset "union_reps! and substitution queueing" begin
    manager = PC.PropagationManager([10, 20])
    @test PC.union_reps!(manager, [2, 3]) == 2
    @test manager.parent_pos[2] == 2
    @test manager.parent_pos[3] == 2
    @test PC.pop_substitution!(manager) == (20, 10, true)
    @test PC.pop_substitution!(manager) === nothing

    manager2 = PC.PropagationManager([1])
    @test PC.union_reps!(manager2, [1, 2]) == 1
    @test PC.pop_substitution!(manager2) === nothing
end

@testset "substitute_scc_by_new_var! swaps SCC representatives to a fresh variable" begin
    manager = prepare_substitution_manager()

    pos_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, false)])]
    neg_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, true)])]
    pos3_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(3, false)])]
    neg3_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(3, true)])]
    pos_component = [manager.lit_to_pos[PC.VarLit(1, false)], manager.lit_to_pos[PC.VarLit(2, false)]]
    neg_component = [manager.lit_to_pos[PC.VarLit(1, true)], manager.lit_to_pos[PC.VarLit(2, true)]]
    split_positions = vcat(pos_component, neg_component)

    @test has_edge(manager.graph, pos_scc, neg3_scc)
    @test has_edge(manager.graph, pos3_scc, pos_scc)
    @test has_edge(manager.graph, pos3_scc, neg_scc)
    @test has_edge(manager.graph, neg_scc, neg3_scc)

    old_ne = ne(manager.graph)
    old_nv = nv(manager.graph)
    old_nreps = manager.nreps

    @test PC.substitute_scc_by_new_var!(manager, 1, 10) === manager
    @test manager.nvars == 4
    @test manager.nreps == old_nreps + length(split_positions)
    @test nv(manager.graph) == old_nv + length(split_positions)
    @test ne(manager.graph) == old_ne
    @test manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(10, false)]] == pos_scc
    @test manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(10, true)]] == neg_scc
    @test manager.scc_pos_to_rep_pos[pos_scc] == manager.lit_to_pos[PC.VarLit(10, false)]
    @test manager.scc_pos_to_rep_pos[neg_scc] == manager.lit_to_pos[PC.VarLit(10, true)]
    @test PC.fixed_value(manager, 10) === nothing
    @test isempty(drain_substitutions!(manager))
    @test isempty(drain_fixings!(manager))

    @test has_edge(manager.graph, pos_scc, neg3_scc)
    @test has_edge(manager.graph, pos3_scc, pos_scc)
    @test has_edge(manager.graph, pos3_scc, neg_scc)
    @test has_edge(manager.graph, neg_scc, neg3_scc)

    isolated_sccs = Int[]
    for pos in split_positions
        @test manager.parent_pos[pos] == pos
        scc = manager.rep_pos_to_scc_pos[pos]
        push!(isolated_sccs, scc)
        @test scc > old_nreps
        @test manager.scc_pos_to_rep_pos[scc] == pos
        @test isempty(outneighbors(manager.graph, scc))
        @test isempty(inneighbors(manager.graph, scc))
    end
    @test length(unique(isolated_sccs)) == length(split_positions)

    @test PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, false)]) == manager.lit_to_pos[PC.VarLit(1, false)]
    @test PC.repr!(manager, manager.lit_to_pos[PC.VarLit(2, false)]) == manager.lit_to_pos[PC.VarLit(2, false)]
    @test PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, true)]) == manager.lit_to_pos[PC.VarLit(1, true)]
    @test PC.repr!(manager, manager.lit_to_pos[PC.VarLit(2, true)]) == manager.lit_to_pos[PC.VarLit(2, true)]
end

@testset "substitute_scc_by_new_var! enforces preconditions" begin
    present_vid_manager = prepare_substitution_manager()
    @test_throws AssertionError PC.substitute_scc_by_new_var!(present_vid_manager, 1, 3)

    labeled_manager = prepare_substitution_manager()
    pos_scc = labeled_manager.rep_pos_to_scc_pos[PC.repr!(labeled_manager, labeled_manager.lit_to_pos[PC.VarLit(1, false)])]
    PC.set_lit_label!(labeled_manager, pos_scc, PC.TRUE)
    @test_throws AssertionError PC.substitute_scc_by_new_var!(labeled_manager, 1, 10)

    singleton_manager = PC.PropagationManager([1])
    @test_throws AssertionError PC.substitute_scc_by_new_var!(singleton_manager, 1, 10)

    collapsed_manager = PC.PropagationManager([1])
    PC.union_reps!(collapsed_manager, [1, 2])
    @test_throws AssertionError PC.substitute_scc_by_new_var!(collapsed_manager, 1, 10)

    queued_manager = prepare_substitution_manager()
    PC.enqueue_substitution!(queued_manager, 99)
    @test_throws AssertionError PC.substitute_scc_by_new_var!(queued_manager, 1, 10)
end

@testset "finalize_phase! strips labeled SCCs and preserves unlabeled structure" begin
    manager = PC.PropagationManager([1, 2, 3, 4])
    PC.add_equivalence!(manager, PC.VarLit(1, false), PC.VarLit(2, false))
    PC.update_sccs!(manager)
    @test Set(drain_substitutions!(manager)) == Set([(2, 1, false)])

    PC.add_implication!(manager, PC.VarLit(3, false), PC.VarLit(4, false))
    PC.add_implication!(manager, PC.VarLit(3, false), PC.VarLit(1, false))
    PC.add_implication!(manager, PC.VarLit(1, false), PC.VarLit(4, false))
    pos3_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(3, false)])]
    pos4_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(4, false)])]
    neg4_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(4, true)])]
    neg3_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(3, true)])]

    pos12_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, false)])
    neg12_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, true)])
    pos12_scc = manager.rep_pos_to_scc_pos[pos12_rep]
    neg12_scc = manager.rep_pos_to_scc_pos[neg12_rep]
    @test has_edge(manager.graph, pos3_scc, pos4_scc)
    @test has_edge(manager.graph, neg4_scc, neg3_scc)
    @test has_edge(manager.graph, pos3_scc, pos12_scc)
    @test has_edge(manager.graph, pos12_scc, pos4_scc)

    PC.fix_var!(manager, 1, true)
    @test Set(drain_fixings!(manager)) == Set([(1, true)])

    old_nreps = manager.nreps
    @test PC.finalize_phase!(manager) === manager
    @test isempty(drain_fixings!(manager))
    @test isempty(drain_substitutions!(manager))
    @test manager.nreps == old_nreps + 2
    @test nv(manager.graph) == manager.nreps
    @test !has_edge(manager.graph, pos3_scc, pos12_scc)
    @test !has_edge(manager.graph, pos12_scc, pos4_scc)
    @test has_edge(manager.graph, pos3_scc, pos4_scc)
    @test has_edge(manager.graph, neg4_scc, neg3_scc)

    @test manager.lit_labels[pos12_scc] == PC.UNDEF
    @test manager.lit_labels[neg12_scc] == PC.UNDEF

    pos_component = [manager.lit_to_pos[PC.VarLit(1, false)], manager.lit_to_pos[PC.VarLit(2, false)]]
    neg_component = [manager.lit_to_pos[PC.VarLit(1, true)], manager.lit_to_pos[PC.VarLit(2, true)]]
    for pos in vcat(pos_component, neg_component)
        @test manager.parent_pos[pos] == pos
        scc = manager.rep_pos_to_scc_pos[pos]
        @test manager.scc_pos_to_rep_pos[scc] == pos
        @test manager.lit_labels[scc] == PC.UNDEF
    end
    @test length(unique(manager.rep_pos_to_scc_pos[pos] for pos in pos_component)) == 2
    @test length(unique(manager.rep_pos_to_scc_pos[pos] for pos in neg_component)) == 2
end

@testset "finalize_phase! and ensure_literals! keep the manager reusable across phases" begin
    manager = PC.PropagationManager([1])
    PC.fix_var!(manager, 1, true)
    @test Set(drain_fixings!(manager)) == Set([(1, true)])
    PC.finalize_phase!(manager)

    @test PC.fixed_value(manager, 1) === nothing
    @test PC.ensure_literals!(manager, [1, 2]) === manager
    PC.add_equivalence!(manager, PC.VarLit(1, false), PC.VarLit(2, false))
    PC.update!(manager)

    @test Set(drain_substitutions!(manager)) == Set([(2, 1, false)])
    @test PC.fixed_value(manager, 1) === nothing
    @test PC.fixed_value(manager, 2) === nothing
end

@testset "Edge and implication builders" begin
    manager = PC.PropagationManager([1, 2])

    @test PC.add_edge!(manager, PC.VarLit(1, false), PC.VarLit(2, true))
    @test has_edge(manager.graph, 1, 4)
    @test !PC.add_edge!(manager, PC.VarLit(1, false), PC.VarLit(2, true))

    manager2 = PC.PropagationManager([1, 2])
    PC.union_reps!(manager2, [1, 3])
    @test !PC.add_edge!(manager2, PC.VarLit(1, false), PC.VarLit(2, false))
    @test ne(manager2.graph) == 0

    manager3 = PC.PropagationManager([1, 2])
    PC.add_implication!(manager3, PC.VarLit(1, false), PC.VarLit(2, true))
    @test has_edge(manager3.graph, 1, 4)
    @test has_edge(manager3.graph, 3, 2)
    @test ne(manager3.graph) == 2

    manager4 = PC.PropagationManager([1, 2])
    PC.add_equivalence!(manager4, PC.VarLit(1, false), PC.VarLit(2, true))
    @test has_edge(manager4.graph, 1, 4)
    @test has_edge(manager4.graph, 3, 2)
    @test has_edge(manager4.graph, 4, 1)
    @test has_edge(manager4.graph, 2, 3)
    @test ne(manager4.graph) == 4
end

@testset "fix_var! labels and queues exactly once" begin
    manager = PC.PropagationManager([5])

    PC.fix_var!(manager, 5, true)
    @test has_edge(manager.graph, 2, 1)
    @test !has_edge(manager.graph, 1, 2)
    @test manager.lit_labels[1] == PC.TRUE
    @test manager.lit_labels[2] == PC.FALSE
    @test PC.pop_fixing!(manager) == (5, true)
    @test PC.pop_fixing!(manager) === nothing

    PC.fix_var!(manager, 5, true)
    @test PC.pop_fixing!(manager) === nothing

    manager2 = PC.PropagationManager([6])
    PC.fix_var!(manager2, 6, false)
    @test has_edge(manager2.graph, 1, 2)
    @test !has_edge(manager2.graph, 2, 1)
    @test manager2.lit_labels[1] == PC.FALSE
    @test manager2.lit_labels[2] == PC.TRUE
    @test PC.pop_fixing!(manager2) == (6, false)
end

@testset "propagation manager keeps the p11/p1 parity relation negated" begin
    pos_to_var = [1, 7, 11, 21]
    conj = falses(4, 4)
    conj[1, 4] = true
    conj[4, 1] = true
    cons = [
        PC.XorConstraint(BitVector([1, 0, 0, 1]), false),
        PC.XorConstraint(falses(4), conj, false),
        PC.XorConstraint(BitVector([1, 1, 0, 0]), true),
        PC.XorConstraint(BitVector([0, 1, 0, 0]), true),
        PC.XorConstraint(BitVector([1, 0, 1, 0]), true),
    ]
    manager = PC.PropagationManager(pos_to_var)

    for con in cons
        PC.register_implications!(manager, con, pos_to_var)
    end
    PC.update!(manager)

    subs = Set(drain_substitutions!(manager))
    @test (11, 1, true) in subs
    @test !((11, 1, false) in subs)
end

@testset "register_implications! adds two-term pure-conjunction cross implications" begin
    pos_to_var = [11, 12, 21, 22]
    con = PC.XorConstraint(
        falses(4),
        propagation_edge_matrix(4, [(1, 2), (3, 4)]),
        true,
    )
    manager = PC.PropagationManager(pos_to_var)

    PC.register_implications!(manager, con, pos_to_var)

    PC.fix_var!(manager, 11, false)
    @test PC.pop_fixing!(manager) == (11, false)
    @test PC.update!(manager) === manager
    @test PC.fixed_value(manager, 21) == true
    @test PC.fixed_value(manager, 22) == true
    @test Set(drain_fixings!(manager)) == Set([(21, true), (22, true)])

    reverse_manager = PC.PropagationManager(pos_to_var)
    reverse_con = PC.XorConstraint(
        falses(4),
        propagation_edge_matrix(4, [(1, 2), (3, 4)]),
        true,
    )

    PC.register_implications!(reverse_manager, reverse_con, pos_to_var)

    PC.fix_var!(reverse_manager, 22, false)
    @test PC.pop_fixing!(reverse_manager) == (22, false)
    @test PC.update!(reverse_manager) === reverse_manager
    @test PC.fixed_value(reverse_manager, 11) == true
    @test PC.fixed_value(reverse_manager, 12) == true
    @test Set(drain_fixings!(reverse_manager)) == Set([(11, true), (12, true)])
end

@testset "register_implications! adds triangle implications for rhs=true" begin
    pos_to_var = [1, 2, 3]
    con = PC.XorConstraint(
        falses(3),
        propagation_edge_matrix(3, [(1, 2), (2, 3), (1, 3)]),
        true,
    )
    manager = PC.PropagationManager(pos_to_var)

    PC.register_implications!(manager, con, pos_to_var)

    PC.fix_var!(manager, 1, false)
    @test PC.pop_fixing!(manager) == (1, false)
    @test PC.update!(manager) === manager
    @test PC.fixed_value(manager, 2) == true
    @test PC.fixed_value(manager, 3) == true
    @test Set(drain_fixings!(manager)) == Set([(2, true), (3, true)])
end

@testset "register_implications! adds triangle implications for rhs=false" begin
    pos_to_var = [1, 2, 3]
    con = PC.XorConstraint(
        falses(3),
        propagation_edge_matrix(3, [(1, 2), (2, 3), (1, 3)]),
        false,
    )
    manager = PC.PropagationManager(pos_to_var)

    PC.register_implications!(manager, con, pos_to_var)

    PC.fix_var!(manager, 1, true)
    @test PC.pop_fixing!(manager) == (1, true)
    @test PC.update!(manager) === manager
    @test PC.fixed_value(manager, 2) == false
    @test PC.fixed_value(manager, 3) == false
    @test Set(drain_fixings!(manager)) == Set([(2, false), (3, false)])
end

@testset "register_implications! skips non-triangle nnz_conj==3 rows" begin
    pos_to_var = [1, 2, 3, 4]
    con = PC.XorConstraint(
        falses(4),
        propagation_edge_matrix(4, [(1, 2), (2, 3), (3, 4)]),
        true,
    )
    manager = PC.PropagationManager(pos_to_var)

    PC.register_implications!(manager, con, pos_to_var)

    PC.fix_var!(manager, 1, false)
    @test PC.pop_fixing!(manager) == (1, false)
    @test PC.update!(manager) === manager
    @test PC.fixed_value(manager, 2) === nothing
    @test PC.fixed_value(manager, 3) === nothing
    @test PC.fixed_value(manager, 4) === nothing
    @test isempty(drain_fixings!(manager))
end

@testset "register_implications! rejects nnz_conj==3 rows whose candidate triangle columns do not match" begin
    pos_to_var = [1, 2, 3, 4]
    con = PC.XorConstraint(
        falses(4),
        propagation_edge_matrix(4, [(1, 2), (1, 3), (2, 4)]),
        true,
    )
    manager = PC.PropagationManager(pos_to_var)

    PC.register_implications!(manager, con, pos_to_var)

    PC.fix_var!(manager, 1, false)
    @test PC.pop_fixing!(manager) == (1, false)
    @test PC.update!(manager) === manager
    @test PC.fixed_value(manager, 2) === nothing
    @test PC.fixed_value(manager, 3) === nothing
    @test PC.fixed_value(manager, 4) === nothing
    @test isempty(drain_fixings!(manager))
end

@testset "update_sccs! condenses graph and carries labels" begin
    manager = PC.PropagationManager([10, 20, 30])
    PC.fix_var!(manager, 10, true)
    @test PC.pop_fixing!(manager) == (10, true)

    PC.add_equivalence!(manager, PC.VarLit(10, false), PC.VarLit(20, false))
    @test ne(manager.graph) == 5

    @test PC.update_sccs!(manager) === manager
    @test manager.nreps == 4
    @test nv(manager.graph) == 4
    @test ne(manager.graph) == 1

    pos10_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(10, false)])
    pos20_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(20, false)])
    neg10_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(10, true)])
    neg20_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(20, true)])

    @test pos10_rep == pos20_rep
    @test neg10_rep == neg20_rep
    @test manager.lit_labels[manager.rep_pos_to_scc_pos[pos10_rep]] == PC.TRUE
    @test manager.lit_labels[manager.rep_pos_to_scc_pos[neg10_rep]] == PC.FALSE
    @test PC.pop_substitution!(manager) == (20, 10, false)
    @test PC.pop_substitution!(manager) === nothing
end

@testset "propagate_labels! pushes labels forward and backward" begin
    manager = PC.PropagationManager([1, 2, 3])
    PC.fix_var!(manager, 1, true)
    @test PC.pop_fixing!(manager) == (1, true)

    PC.add_implication!(manager, PC.VarLit(1, false), PC.VarLit(2, false))
    PC.add_implication!(manager, PC.VarLit(2, false), PC.VarLit(3, true))

    top_order = topological_sort(manager.graph)
    @test PC.propagate_labels!(manager, top_order) === manager

    @test manager.lit_labels[manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(2, false)]]] == PC.TRUE
    @test manager.lit_labels[manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(2, true)]]] == PC.FALSE
    @test manager.lit_labels[manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(3, false)]]] == PC.FALSE
    @test manager.lit_labels[manager.rep_pos_to_scc_pos[manager.lit_to_pos[PC.VarLit(3, true)]]] == PC.TRUE

    fixings = Set(drain_fixings!(manager))
    @test fixings == Set([(2, true), (3, false)])
end

@testset "label_reachables! finds forced literals" begin
    manager = PC.PropagationManager([1, 2])
    PC.add_edge!(manager, PC.VarLit(1, true), PC.VarLit(1, false))
    PC.add_edge!(manager, PC.VarLit(2, false), PC.VarLit(2, true))

    top_order = topological_sort(manager.graph)
    @test PC.label_reachables!(manager, top_order) === manager

    @test manager.lit_labels[1] == PC.TRUE
    @test manager.lit_labels[2] == PC.FALSE
    @test manager.lit_labels[3] == PC.FALSE
    @test manager.lit_labels[4] == PC.TRUE

    fixings = Set(drain_fixings!(manager))
    @test fixings == Set([(1, true), (2, false)])

    PC.label_reachables!(manager, top_order)
    @test PC.pop_fixing!(manager) === nothing
end

@testset "Queue pop functions are FIFO" begin
    manager = PC.PropagationManager([1, 2])
    PC.enqueue_fixing!(manager, 1, true)
    PC.enqueue_fixing!(manager, 2, false)
    @test PC.pop_fixing!(manager) == (1, true)
    @test PC.pop_fixing!(manager) == (2, false)
    @test PC.pop_fixing!(manager) === nothing

    PC.enqueue_substitution!(manager, 1)
    PC.enqueue_substitution!(manager, 2)
    @test PC.pop_substitution!(manager) == (1, 1, false)
    @test PC.pop_substitution!(manager) == (2, 2, false)
    @test PC.pop_substitution!(manager) === nothing
end

@testset "update! integrates SCC updates, reachability, and propagation" begin
    manager = PC.PropagationManager([1, 2, 3])
    PC.fix_var!(manager, 1, true)
    PC.add_equivalence!(manager, PC.VarLit(1, false), PC.VarLit(2, false))
    PC.add_implication!(manager, PC.VarLit(2, false), PC.VarLit(3, false))

    @test PC.update!(manager) === manager
    @test manager.nreps == 4

    pos1_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, false)])
    pos2_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(2, false)])
    neg1_rep = PC.repr!(manager, manager.lit_to_pos[PC.VarLit(1, true)])
    pos3_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(3, false)])]
    neg3_scc = manager.rep_pos_to_scc_pos[PC.repr!(manager, manager.lit_to_pos[PC.VarLit(3, true)])]
    pos12_scc = manager.rep_pos_to_scc_pos[pos1_rep]
    neg12_scc = manager.rep_pos_to_scc_pos[neg1_rep]

    @test pos1_rep == pos2_rep
    @test manager.lit_labels[pos12_scc] == PC.TRUE
    @test manager.lit_labels[neg12_scc] == PC.FALSE
    @test manager.lit_labels[pos3_scc] == PC.TRUE
    @test manager.lit_labels[neg3_scc] == PC.FALSE

    @test Set(drain_fixings!(manager)) == Set([(1, true), (3, true)])
    @test PC.pop_substitution!(manager) == (2, 1, false)
    @test PC.pop_substitution!(manager) === nothing
end
