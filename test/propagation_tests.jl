using Test
using Graphs: has_edge, ne, nv, topological_sort
import QIPresolve.PresolvingCore as PC

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
    @test has_edge(manager.graph, 1, 2)
    @test manager.lit_labels[1] == PC.TRUE
    @test manager.lit_labels[2] == PC.FALSE
    @test PC.pop_fixing!(manager) == (5, true)
    @test PC.pop_fixing!(manager) === nothing

    PC.fix_var!(manager, 5, true)
    @test PC.pop_fixing!(manager) === nothing

    manager2 = PC.PropagationManager([6])
    PC.fix_var!(manager2, 6, false)
    @test has_edge(manager2.graph, 2, 1)
    @test manager2.lit_labels[1] == PC.FALSE
    @test manager2.lit_labels[2] == PC.TRUE
    @test PC.pop_fixing!(manager2) == (6, false)
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
    @test manager.lit_labels[pos12_scc] == PC.FALSE
    @test manager.lit_labels[neg12_scc] == PC.TRUE
    @test manager.lit_labels[pos3_scc] == PC.UNDEF
    @test manager.lit_labels[neg3_scc] == PC.UNDEF

    @test Set(drain_fixings!(manager)) == Set([(1, true), (1, false)])
    @test PC.pop_substitution!(manager) == (2, 1, false)
    @test PC.pop_substitution!(manager) === nothing
end
