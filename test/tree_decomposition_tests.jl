using Test
import Graphs
using Graphs: SimpleGraph
import QIPresolve.PresolvingCore as PC

function td_component(var_ids::Vector{PC.VarId}, graph_edges::Vector{Tuple{Int, Int}})
    graph = SimpleGraph{Int}(length(var_ids))
    for (source, destination) in graph_edges
        Graphs.add_edge!(graph, source, destination)
    end

    return PC.NonSingleton(
        graph,
        Dict(var_id => position for (position, var_id) in enumerate(var_ids)),
        var_ids,
        Dict{PC.VarId, Int}(),
        Dict{Tuple{PC.VarId, PC.VarId}, Int}(),
    )
end

function td_tree(td::PC.TreeDecomposition)
    tree = SimpleGraph{Int}(length(td.bags))
    roots = Int[]
    for node in eachindex(td.parents)
        parent = td.parents[node]
        if parent == 0
            push!(roots, node)
        else
            Graphs.add_edge!(tree, node, parent)
        end
    end

    @test length(roots) == 1
    @test Graphs.ne(tree) == Graphs.nv(tree) - 1
    @test Graphs.is_connected(tree)
    return tree
end

function td_tree(nice::PC.NiceTreeDecomposition)
    tree = SimpleGraph{Int}(length(nice.bags))
    parent_count = zeros(Int, length(nice.bags))
    for node in eachindex(nice.children)
        for child in nice.children[node]
            Graphs.add_edge!(tree, node, child)
            parent_count[child] += 1
        end
    end

    @test parent_count[nice.root] == 0
    @test all(
        parent_count[node] == 1
        for node in eachindex(parent_count)
        if node != nice.root
    )
    @test Graphs.ne(tree) == Graphs.nv(tree) - 1
    @test Graphs.is_connected(tree)
    return tree
end

function assert_tree_decomposition(graph::SimpleGraph, bags::Vector{BitSet}, tree::SimpleGraph)
    covered_vertices = BitSet()
    for bag in bags
        union!(covered_vertices, bag)
    end
    @test covered_vertices == BitSet(1:Graphs.nv(graph))

    for edge in Graphs.edges(graph)
        source = Graphs.src(edge)
        destination = Graphs.dst(edge)
        @test any(source in bag && destination in bag for bag in bags)
    end

    for vertex in Graphs.vertices(graph)
        containing_nodes = BitSet(
            node for node in eachindex(bags) if vertex in bags[node]
        )
        @test !isempty(containing_nodes)

        seen = BitSet()
        stack = Int[first(containing_nodes)]
        while !isempty(stack)
            node = pop!(stack)
            node in seen && continue
            push!(seen, node)
            for neighbor in Graphs.neighbors(tree, node)
                neighbor in containing_nodes && push!(stack, neighbor)
            end
        end
        @test seen == containing_nodes
    end

    return nothing
end

td_width(bags::Vector{BitSet}) = maximum(length, bags) - 1

function assert_nice_node_types(
        nice::PC.NiceTreeDecomposition,
        component::PC.NonSingleton,
    )
    for node in eachindex(nice.bags)
        bag = nice.bags[node]
        children = nice.children[node]
        action = nice.actions[node]

        if action isa PC.LeafAction
            @test isempty(children)
            @test isempty(bag)
        elseif action isa PC.IntroduceAction
            @test length(children) == 1
            child_bag = nice.bags[only(children)]
            @test length(bag) == length(child_bag) + 1
            @test issubset(child_bag, bag)
            introduced_position = component.var_id_to_pos[action.id]
            @test setdiff(bag, child_bag) == BitSet([introduced_position])
        elseif action isa PC.RemoveAction
            @test length(children) == 1
            child_bag = nice.bags[only(children)]
            @test length(child_bag) == length(bag) + 1
            @test issubset(bag, child_bag)
            removed_position = component.var_id_to_pos[action.id]
            @test setdiff(child_bag, bag) == BitSet([removed_position])
        elseif action isa PC.JoinAction
            @test length(children) == 2
            @test bag == nice.bags[children[1]]
            @test bag == nice.bags[children[2]]
        else
            @test false
        end
    end

    @test isempty(nice.bags[nice.root])
    @test all(
        isempty(nice.bags[node])
        for node in eachindex(nice.children)
        if isempty(nice.children[node])
    )

    for node in eachindex(nice.children)
        for child in nice.children[node]
            difference = symdiff(nice.bags[node], nice.bags[child])
            @test length(difference) <= 1
        end
    end

    return nothing
end

function interpret_action_order(
        actions::Vector{PC.ResidueAction},
        expected_var_ids::Vector{PC.VarId},
    )
    states = Set{PC.VarId}[]
    remove_counts = Dict(var_id => 0 for var_id in expected_var_ids)

    for action in actions
        if action isa PC.LeafAction
            push!(states, Set{PC.VarId}())
        elseif action isa PC.IntroduceAction
            @test !isempty(states)
            state = pop!(states)
            @test !(action.id in state)
            push!(state, action.id)
            push!(states, state)
        elseif action isa PC.RemoveAction
            @test !isempty(states)
            state = pop!(states)
            @test action.id in state
            delete!(state, action.id)
            remove_counts[action.id] = get(remove_counts, action.id, 0) + 1
            push!(states, state)
        elseif action isa PC.JoinAction
            @test length(states) >= 2
            right_state = pop!(states)
            left_state = pop!(states)
            @test left_state == right_state
            push!(states, left_state)
        else
            @test false
        end
    end

    @test length(states) == 1
    @test isempty(only(states))
    @test remove_counts == Dict(var_id => 1 for var_id in expected_var_ids)
    return nothing
end

function action_signature(action::PC.ResidueAction)
    if action isa PC.IntroduceAction
        return (:introduce, action.id)
    elseif action isa PC.RemoveAction
        return (:remove, action.id)
    elseif action isa PC.LeafAction
        return (:leaf, nothing)
    end
    return (:join, nothing)
end

@testset "minimum-degree tree decomposition" begin
    path = td_component(
        PC.VarId[10, 20, 30, 50],
        Tuple{Int, Int}[(1, 2), (2, 3), (3, 4)],
    )
    cycle = td_component(
        PC.VarId[2, 7, 11, 20],
        Tuple{Int, Int}[(1, 2), (2, 3), (3, 4), (4, 1)],
    )
    star = td_component(
        PC.VarId[5, 10, 15, 20, 25],
        Tuple{Int, Int}[(3, 1), (3, 2), (3, 4), (3, 5)],
    )

    for component in (path, cycle, star)
        td = PC.minimum_degree_tree_decomposition(component)
        tree = td_tree(td)
        assert_tree_decomposition(component.graph, td.bags, tree)
        @test td.parents[end] == 0
    end

    cycle_td = PC.minimum_degree_tree_decomposition(cycle)
    @test cycle_td.elimination_order == [1, 2, 3, 4]
    @test td_width(cycle_td.bags) == 2
    @test Graphs.ne(cycle.graph) == 4
    @test !Graphs.has_edge(cycle.graph, 2, 4)

    star_td = PC.minimum_degree_tree_decomposition(star)
    @test star_td.elimination_order == [1, 2, 4, 3, 5]
end

@testset "nice tree decomposition preserves validity and width" begin
    components = [
        td_component(
            PC.VarId[10, 20, 30, 50],
            Tuple{Int, Int}[(1, 2), (2, 3), (3, 4)],
        ),
        td_component(
            PC.VarId[2, 7, 11, 20],
            Tuple{Int, Int}[(1, 2), (2, 3), (3, 4), (4, 1)],
        ),
        td_component(
            PC.VarId[5, 10, 15, 20, 25],
            Tuple{Int, Int}[(3, 1), (3, 2), (3, 4), (3, 5)],
        ),
    ]

    for component in components
        td = PC.minimum_degree_tree_decomposition(component)
        nice = PC.make_nice_tree_decomposition(td, component.pos_to_var_id)

        nice_tree = td_tree(nice)
        assert_tree_decomposition(component.graph, nice.bags, nice_tree)
        assert_nice_node_types(nice, component)
        @test td_width(nice.bags) == td_width(td.bags)

        actions = PC.action_order(component)
        @test length(actions) == length(nice.bags)
        interpret_action_order(actions, component.pos_to_var_id)
    end

    star = components[3]
    star_actions = PC.action_order(star)
    @test any(action isa PC.JoinAction for action in star_actions)
    @test action_signature.(star_actions) ==
          action_signature.(PC.action_order(star))
end

@testset "action_order validates NonSingleton invariants" begin
    singleton = td_component(PC.VarId[1], Tuple{Int, Int}[])
    disconnected = td_component(PC.VarId[1, 2], Tuple{Int, Int}[])

    @test_throws ArgumentError PC.action_order(singleton)
    @test_throws ArgumentError PC.action_order(disconnected)
end
