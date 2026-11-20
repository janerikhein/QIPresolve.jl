using DataStructures: PriorityQueue, dequeue!

"""
    ResidueAction

Abstract supertype for bottom-up operations on a nice tree decomposition.
"""
abstract type ResidueAction end

"""
    LeafAction()

Create the empty state associated with a leaf of a nice tree decomposition.
"""
struct LeafAction <: ResidueAction end

"""
    IntroduceAction(id)

Introduce variable `id` into the state produced by the node's child.
"""
struct IntroduceAction <: ResidueAction
    id::VarId
end

"""
    RemoveAction(id)

Remove variable `id` from the state produced by the node's child.
"""
struct RemoveAction <: ResidueAction
    id::VarId
end

"""
    JoinAction()

Join the two states produced by the children of a join node.
"""
struct JoinAction <: ResidueAction end

"""
    TreeDecomposition

Tree decomposition whose nodes are indexed by the entries of `bags`.
`parents[root]` is zero, and all other entries contain the parent node index.
Bag elements and `elimination_order` are positions in the source graph.
"""
struct TreeDecomposition
    bags::Vector{BitSet}
    parents::Vector{Int}
    elimination_order::Vector{Int}
end

"""
    NiceTreeDecomposition

Rooted nice tree decomposition. Each node has a bag, an ordered list of zero,
one, or two children, and the corresponding bottom-up action.
"""
struct NiceTreeDecomposition
    bags::Vector{BitSet}
    children::Vector{Vector{Int}}
    root::Int
    actions::Vector{ResidueAction}
end

tree_decomposition_width(td::TreeDecomposition) = maximum(length, td.bags; init = 0) - 1

function _validate_non_singleton(component::NonSingleton)
    nvertices = Graphs.nv(component.graph)
    nvertices >= 2 ||
        throw(ArgumentError("NonSingleton graph must contain at least two vertices"))
    Graphs.is_connected(component.graph) ||
        throw(ArgumentError("NonSingleton graph must be connected"))
    length(component.pos_to_var_id) == nvertices ||
        throw(ArgumentError("pos_to_var_id must contain one entry per graph vertex"))

    return nothing
end

"""
    minimum_degree_tree_decomposition(component)

Construct a deterministic tree decomposition of `component.graph` using the
minimum-degree elimination heuristic. Degree ties are broken by variable id.
"""
function minimum_degree_tree_decomposition(component::NonSingleton)
    _validate_non_singleton(component)

    nvertices = Graphs.nv(component.graph)
    adjacency = [
        BitSet(Graphs.neighbors(component.graph, vertex))
        for vertex in 1:nvertices
    ]

    queue = PriorityQueue{Int, Tuple{Int, VarId}}()
    for vertex in 1:nvertices
        queue[vertex] = (length(adjacency[vertex]), component.pos_to_var_id[vertex])
    end

    elimination_order = Int[]
    bags = BitSet[]
    remaining_neighbors = BitSet[]
    sizehint!(elimination_order, nvertices)
    sizehint!(bags, nvertices)
    sizehint!(remaining_neighbors, nvertices)

    while !isempty(queue)
        vertex = dequeue!(queue)
        neighbors = collect(adjacency[vertex])

        bag = copy(adjacency[vertex])
        push!(bag, vertex)
        push!(bags, bag)
        push!(remaining_neighbors, copy(adjacency[vertex]))
        push!(elimination_order, vertex)

        for first_index in eachindex(neighbors)
            first_neighbor = neighbors[first_index]
            for second_index in (first_index + 1):length(neighbors)
                second_neighbor = neighbors[second_index]
                second_neighbor in adjacency[first_neighbor] && continue
                push!(adjacency[first_neighbor], second_neighbor)
                push!(adjacency[second_neighbor], first_neighbor)
            end
        end

        for neighbor in neighbors
            delete!(adjacency[neighbor], vertex)
            queue[neighbor] = (
                length(adjacency[neighbor]),
                component.pos_to_var_id[neighbor],
            )
        end
        empty!(adjacency[vertex])
    end

    elimination_position = zeros(Int, nvertices)
    for (position, vertex) in enumerate(elimination_order)
        elimination_position[vertex] = position
    end

    parents = zeros(Int, nvertices)
    for node in 1:(nvertices - 1)
        neighbors = remaining_neighbors[node]
        isempty(neighbors) &&
            error("connected elimination graph has no parent candidate for node $node")

        parent_vertex = first(neighbors)
        parent_position = elimination_position[parent_vertex]
        for candidate in Iterators.drop(neighbors, 1)
            candidate_position = elimination_position[candidate]
            if candidate_position < parent_position
                parent_vertex = candidate
                parent_position = candidate_position
            end
        end
        parents[node] = parent_position
    end

    return TreeDecomposition(bags, parents, elimination_order)
end

function _postorder(children::Vector{Vector{Int}}, root::Int)
    order = Int[]
    stack = Tuple{Int, Bool}[(root, false)]

    while !isempty(stack)
        node, expanded = pop!(stack)
        if expanded
            push!(order, node)
            continue
        end

        push!(stack, (node, true))
        for child in Iterators.reverse(children[node])
            push!(stack, (child, false))
        end
    end

    return order
end

function _reachable_preorder(children::Vector{Vector{Int}}, root::Int)
    order = Int[]
    stack = Int[root]

    while !isempty(stack)
        node = pop!(stack)
        push!(order, node)
        for child in Iterators.reverse(children[node])
            push!(stack, child)
        end
    end

    return order
end

function _infer_action(
        bag::BitSet,
        children::Vector{Int},
        bags::Vector{BitSet},
        pos_to_var_id::Vector{VarId},
    )::ResidueAction
    if isempty(children)
        isempty(bag) || error("leaf nodes must have empty bags")
        return LeafAction()
    end

    if length(children) == 2
        bag == bags[children[1]] == bags[children[2]] ||
            error("join node and child bags must be identical")
        return JoinAction()
    end

    length(children) == 1 || error("nice tree nodes may have at most two children")
    child_bag = bags[only(children)]

    if length(bag) == length(child_bag) + 1 && issubset(child_bag, bag)
        position = only(setdiff(bag, child_bag))
        return IntroduceAction(pos_to_var_id[position])
    end

    if length(child_bag) == length(bag) + 1 && issubset(bag, child_bag)
        position = only(setdiff(child_bag, bag))
        return RemoveAction(pos_to_var_id[position])
    end

    error("unary nice tree nodes must differ from their child by one vertex")
end

"""
    make_nice_tree_decomposition(td, pos_to_var_id)

Transform `td` into a rooted nice tree decomposition with empty root and leaves.
The conversion contracts equal adjacent bags and never increases the width.
"""
function make_nice_tree_decomposition(
        td::TreeDecomposition,
        pos_to_var_id::Vector{VarId},
    )
    nnodes = length(td.bags)
    nnodes > 0 || throw(ArgumentError("tree decomposition must not be empty"))
    length(td.parents) == nnodes ||
        throw(ArgumentError("tree decomposition must have one parent entry per bag"))

    original_children = [Int[] for _ in 1:nnodes]
    original_root = 0
    for node in 1:nnodes
        parent = td.parents[node]
        if parent == 0
            original_root == 0 ||
                throw(ArgumentError("tree decomposition must have exactly one root"))
            original_root = node
        elseif 1 <= parent <= nnodes
            push!(original_children[parent], node)
        else
            throw(ArgumentError("invalid parent index $parent for tree node $node"))
        end
    end
    original_root != 0 ||
        throw(ArgumentError("tree decomposition must have exactly one root"))

    bags = BitSet[]
    children = Vector{Int}[]

    function new_node(bag::BitSet)
        push!(bags, copy(bag))
        push!(children, Int[])
        return length(bags)
    end

    function add_child(parent::Int, child::Int)
        push!(children[parent], child)
        return nothing
    end

    function connect_with_transition(upper::Int, lower::Int)
        upper_bag = bags[upper]
        lower_bag = bags[lower]

        if upper_bag == lower_bag
            children[upper] = copy(children[lower])
            return nothing
        end

        removed = collect(setdiff(upper_bag, lower_bag))
        added = collect(setdiff(lower_bag, upper_bag))
        sort!(removed; by = position -> pos_to_var_id[position])
        sort!(added; by = position -> pos_to_var_id[position])

        transition_bags = BitSet[]
        current_bag = copy(upper_bag)
        for position in removed
            delete!(current_bag, position)
            push!(transition_bags, copy(current_bag))
        end
        for position in added
            push!(current_bag, position)
            push!(transition_bags, copy(current_bag))
        end

        current = upper
        for index in 1:(length(transition_bags) - 1)
            intermediate = new_node(transition_bags[index])
            add_child(current, intermediate)
            current = intermediate
        end
        add_child(current, lower)
        return nothing
    end

    function create_join_ports(root::Int, count::Int)
        count >= 1 || throw(ArgumentError("join tree requires at least one port"))

        ports = Int[]
        tasks = Tuple{Int, Int}[(root, count)]
        while !isempty(tasks)
            node, node_count = pop!(tasks)
            if node_count == 1
                push!(ports, node)
                continue
            end

            left_count = fld(node_count, 2)
            right_count = node_count - left_count
            left = new_node(bags[node])
            right = new_node(bags[node])
            add_child(node, left)
            add_child(node, right)

            push!(tasks, (right, right_count))
            push!(tasks, (left, left_count))
        end

        return ports
    end

    subtree_roots = zeros(Int, nnodes)
    for node in _postorder(original_children, original_root)
        output_root = new_node(td.bags[node])
        child_roots = [subtree_roots[child] for child in original_children[node]]

        if isempty(child_roots)
            empty_leaf = new_node(BitSet())
            connect_with_transition(output_root, empty_leaf)
        elseif length(child_roots) == 1
            connect_with_transition(output_root, only(child_roots))
        else
            ports = create_join_ports(output_root, length(child_roots))
            for index in eachindex(child_roots)
                connect_with_transition(ports[index], child_roots[index])
            end
        end

        subtree_roots[node] = output_root
    end

    empty_root = new_node(BitSet())
    connect_with_transition(empty_root, subtree_roots[original_root])

    reachable = _reachable_preorder(children, empty_root)
    old_to_new = zeros(Int, length(bags))
    for (new, old) in enumerate(reachable)
        old_to_new[old] = new
    end

    compact_bags = [bags[old] for old in reachable]
    compact_children = [
        [old_to_new[child] for child in children[old]]
        for old in reachable
    ]
    compact_root = old_to_new[empty_root]
    actions = ResidueAction[
        _infer_action(
            compact_bags[node],
            compact_children[node],
            compact_bags,
            pos_to_var_id,
        )
        for node in eachindex(compact_bags)
    ]

    return NiceTreeDecomposition(
        compact_bags,
        compact_children,
        compact_root,
        actions,
    )
end

function _action_order(nice::NiceTreeDecomposition)::Vector{ResidueAction}
    return ResidueAction[
        nice.actions[node]
        for node in _postorder(nice.children, nice.root)
    ]
end

"""
    action_order(component) -> Vector{ResidueAction}

Return the deterministic bottom-up action order induced by a nice tree
decomposition of `component.graph`.
"""
function action_order(component::NonSingleton)::Vector{ResidueAction}
    td = minimum_degree_tree_decomposition(component)
    return action_order(component, td)
end

function action_order(component::NonSingleton, td::TreeDecomposition)::Vector{ResidueAction}
    nice = make_nice_tree_decomposition(td, component.pos_to_var_id)
    return _action_order(nice)
end
