using Graphs

const VarId = Int

"""
    SCCNode

Store one strongly connected component in literal space.

`lits` marks positive literals, and `lits_neg` marks negated literals,
both indexed by variable position `1:nvars`.
"""
struct SCCNode
    lits::BitVector
    lits_neg::BitVector
end

"""
    ParityPropagator

Maintain a literal implication graph together with its SCC condensation.

The propagator uses two coupled representations:
1. A union-find structure (`parent_pos`) over literal positions.
2. A condensed implication graph (`scc_graph`) whose nodes correspond to SCCs,
   with payload bitsets stored in `sccs`.

Literal positions use the convention:
- `1:nvars` for positive literals.
- `(nvars+1):(2*nvars)` for negated literals.

"""
mutable struct ParityPropagator
    nvars::Int

    # Union Find mappings
    var_id_to_var_pos::Dict{VarId, Int}
    var_pos_to_var_id::Vector{VarId}
    parent_pos::Vector{Int}

    # Implication graph mappings
    var_pos_to_node_idx::Dict{Int, Int}
    node_idx_to_var_pos::Vector{Int}
    scc_graph::DiGraph
    sccs::Vector{SCCNode}
    requires_update::Bool
end


function _build_parity_propagator(
    nvars::Int,
    var_id_to_var_pos::Dict{VarId, Int},
    var_pos_to_var_id::Vector{VarId},
)
    nlits = 2 * nvars
    parent_pos = collect(1:nlits)

    var_pos_to_node_idx = Dict{Int, Int}(pos => pos for pos in 1:nlits)
    node_idx_to_var_pos = collect(1:nlits)
    scc_graph = DiGraph(nlits)
    sccs = Vector{SCCNode}(undef, nlits)
    for pos in 1:nlits
        lits = falses(nvars)
        lits_neg = falses(nvars)
        if pos <= nvars
            lits[pos] = true
        else
            lits_neg[pos - nvars] = true
        end
        sccs[pos] = SCCNode(lits, lits_neg)
    end

    return ParityPropagator(
        nvars,
        var_id_to_var_pos,
        var_pos_to_var_id,
        parent_pos,
        var_pos_to_node_idx,
        node_idx_to_var_pos,
        scc_graph,
        sccs,
        false,
    )
end

"""
    ParityPropagator(nvars::Int)

Construct a parity propagator with `2*nvars` singleton literal SCCs.

The created propagator has:
- identity variable/literal mappings for `1:nvars`
- an empty directed implication graph with `2*nvars` nodes
- one `SCCNode` per literal
- `requires_update == false`

# Arguments
- `nvars::Int`: Number of variables.

# Throws
- `ArgumentError`: If `nvars < 0`.
"""
function ParityPropagator(nvars::Int)
    nvars >= 0 || throw(ArgumentError("nvars must be non-negative"))

    var_pos_to_var_id = collect(1:nvars)
    var_id_to_var_pos = Dict{VarId, Int}(vid => pos for (pos, vid) in enumerate(var_pos_to_var_id))
    return _build_parity_propagator(nvars, var_id_to_var_pos, var_pos_to_var_id)
end

"""
    ParityPropagator(var_id_to_var_pos::Dict{VarId,Int}, var_pos_to_var_id::Vector{VarId})

Construct a parity propagator with explicit variable-position mappings.

This is useful when literal positions must match an existing parity model
ordering.

# Throws
- `ArgumentError`: If mappings are inconsistent or not a bijection over `1:nvars`.
"""
function ParityPropagator(var_id_to_var_pos::Dict{VarId, Int}, var_pos_to_var_id::Vector{VarId})
    nvars = length(var_pos_to_var_id)
    length(var_id_to_var_pos) == nvars ||
        throw(ArgumentError("var_id_to_var_pos and var_pos_to_var_id must have same length"))

    seen_pos = falses(nvars)
    for (vid, pos) in var_id_to_var_pos
        1 <= pos <= nvars || throw(ArgumentError("invalid position $pos for variable $vid"))
        var_pos_to_var_id[pos] == vid ||
            throw(ArgumentError("mapping mismatch at position $pos: expected $(var_pos_to_var_id[pos]), got $vid"))
        seen_pos[pos] && throw(ArgumentError("duplicate position $pos in var_id_to_var_pos"))
        seen_pos[pos] = true
    end
    all(seen_pos) || throw(ArgumentError("var_id_to_var_pos must cover all positions 1:$nvars"))

    return _build_parity_propagator(nvars, copy(var_id_to_var_pos), copy(var_pos_to_var_id))
end

"""
    lit_to_pos(prop::ParityPropagator, vid::VarId, neg::Bool) -> Int

Map a literal `(vid, neg)` to its internal literal position.

# Arguments
- `prop`: Propagator containing variable-position mapping.
- `vid`: Variable identifier.
- `neg`: `true` for negated literal, `false` for positive literal.

# Returns
- Literal position in `1:(2*prop.nvars)`.
"""
function lit_to_pos(prop::ParityPropagator, vid::VarId, neg::Bool)
    pos = prop.var_id_to_var_pos[vid]
    return neg ? pos + prop.nvars : pos
end

"""
    pos_to_lit(prop::ParityPropagator, pos::Int) -> Tuple{VarId,Bool}

Map an internal literal position back to `(vid, neg)`.

# Arguments
- `prop`: Propagator containing position-variable mapping.
- `pos`: Literal position in `1:(2*prop.nvars)`.

# Returns
- Tuple `(vid, neg)` where `neg` indicates whether the literal is negated.
"""
function pos_to_lit(prop::ParityPropagator, pos::Int)
    neg = pos > prop.nvars
    pos = mod1(pos, prop.nvars)
    vid = prop.var_pos_to_var_id[pos]

    return (vid, neg)
end

"""
    add_edge!(prop::ParityPropagator, vid1::VarId, neg1::Bool, vid2::VarId, neg2::Bool) -> Bool

Add implication edge `(vid1, neg1) -> (vid2, neg2)` to the SCC graph.

Both endpoints are first mapped to their current SCC representatives. If the edge
did not already exist, `requires_update` is set to `true`.

# Returns
- `true` if a new edge was inserted, `false` otherwise.
"""
function add_edge!(prop::ParityPropagator, vid1::VarId, neg1::Bool, vid2::VarId, neg2::Bool)
    pos1 = lit_to_pos(prop, vid1, neg1)
    pos2 = lit_to_pos(prop, vid2, neg2)

    repr1 = get_scc_repr_pos!(prop, pos1)
    repr2 = get_scc_repr_pos!(prop, pos2)
    node1 = prop.var_pos_to_node_idx[repr1]
    node2 = prop.var_pos_to_node_idx[repr2]

    added = Graphs.add_edge!(prop.scc_graph, node1, node2)
    prop.requires_update |= added
    return added
end

"""
    fix_var!(prop::ParityPropagator, vid::VarId, value::Bool) -> Bool

Fix variable `vid` to a Boolean value in the implication graph.

This is encoded as a unit implication:
- `value == true`: add edge `¬vid -> vid`
- `value == false`: add edge `vid -> ¬vid`

The return value matches [`add_edge!`](@ref).
"""
function fix_var!(prop::ParityPropagator, vid::VarId, value::Bool)
    if value
        return add_edge!(prop, vid, true, vid, false)
    else
        return add_edge!(prop, vid, false, vid, true)
    end
end

"""
    update!(prop::ParityPropagator) -> ParityPropagator

Recompute SCCs and condensed implication graph when pending edge updates exist.

If `prop.requires_update == false`, this is a no-op. Otherwise, the function:
1. recomputes SCCs of `prop.scc_graph`
2. rebuilds the condensed DAG between SCCs
3. merges union-find representatives with `union_scc!`
4. rebuilds SCC payloads and node-index mappings

# Returns
- The mutated `prop`.
"""
function update!(prop::ParityPropagator)
    # Nothing changed since last condensation update.
    !prop.requires_update && return prop

    # Snapshot current state.
    nvars = prop.nvars
    old_graph = prop.scc_graph
    old_sccs = prop.sccs
    old_node_idx_to_var_pos = prop.node_idx_to_var_pos

    # Compute SCCs of the current graph and map old nodes -> new SCC index.
    comps = strongly_connected_components(old_graph)
    ncomps = length(comps)
    old_to_new = zeros(Int, nv(old_graph))
    for (new_idx, comp) in enumerate(comps)
        for old_idx in comp
            old_to_new[old_idx] = new_idx
        end
    end

    # Build condensation DAG between SCCs.
    new_graph = DiGraph(ncomps)
    for edge in edges(old_graph)
        src_comp = old_to_new[src(edge)]
        dst_comp = old_to_new[dst(edge)]
        src_comp == dst_comp && continue
        Graphs.add_edge!(new_graph, src_comp, dst_comp)
    end

    # Rebuild SCC payloads and node-index mappings, and merge parents incrementally.
    var_pos_to_node_idx = Dict{Int, Int}()
    node_idx_to_var_pos = Vector{Int}(undef, ncomps)
    sccs = Vector{SCCNode}(undef, ncomps)
    for (new_idx, comp) in enumerate(comps)
        # Merge previous SCC representatives that now belong to one SCC.
        repr_pos = old_node_idx_to_var_pos[first(comp)]
        for old_idx in @view comp[2:end]
            repr_pos = union_scc!(prop, repr_pos, old_node_idx_to_var_pos[old_idx])
        end
        repr_pos = get_scc_repr_pos!(prop, repr_pos)

        # Merge literal bitsets of all old nodes inside this SCC.
        lits = falses(nvars)
        lits_neg = falses(nvars)
        for old_idx in comp
            lits .|= old_sccs[old_idx].lits
            lits_neg .|= old_sccs[old_idx].lits_neg
        end

        # Store SCC data and representative <-> node index mapping.
        sccs[new_idx] = SCCNode(lits, lits_neg)
        node_idx_to_var_pos[new_idx] = repr_pos
        var_pos_to_node_idx[repr_pos] = new_idx
    end

    # Commit rebuilt structures and clear dirty flag.
    prop.var_pos_to_node_idx = var_pos_to_node_idx
    prop.node_idx_to_var_pos = node_idx_to_var_pos
    prop.scc_graph = new_graph
    prop.sccs = sccs
    prop.requires_update = false

    return prop
end

"""
    propagate!(prop::ParityPropagator) -> Tuple{Bool,BitVector,BitVector}

Propagate fixed truth values through SCC implications.

The routine first calls [`update!`](@ref), then applies:
1. SCC consistency check (`x` and `¬x` in same SCC => infeasible).
2. Forced assignments from paths `¬v -> v`.
3. Forward reachability propagation from true literals.
4. Contradiction checks (true and false simultaneously, or edge true -> false).

# Returns
- `(is_feasible, lit_is_true, lit_is_false)` where:
- `is_feasible::Bool` indicates whether a contradiction was found.
- `lit_is_true::BitVector` marks literals forced to true.
- `lit_is_false::BitVector` marks literals forced to false.
- Literal indexing follows `1:nvars` (positive) and `(nvars+1):(2*nvars)` (negated).
"""
function propagate!(prop::ParityPropagator)
    update!(prop)

    nvars = prop.nvars
    nlits = 2 * nvars
    g = prop.scc_graph
    nnodes = nv(g)
    neg_pos(pos::Int) = (pos <= nvars) ? (pos + nvars) : (pos - nvars)
    node_for_pos(pos::Int) = prop.var_pos_to_node_idx[get_scc_repr_pos!(prop, pos)]

    lit_is_true = falses(nlits)
    lit_is_false = falses(nlits)

    # Infeasible if one SCC contains a literal and its negation.
    for scc in prop.sccs
        any(scc.lits .& scc.lits_neg) && return false, lit_is_true, lit_is_false
    end

    # Materialize literals per SCC-node as bitsets over literal positions.
    node_lits = Vector{BitSet}(undef, nnodes)
    for node_idx in 1:nnodes
        lits = BitSet()
        scc = prop.sccs[node_idx]
        for vid in 1:nvars
            scc.lits[vid] && push!(lits, vid)
            scc.lits_neg[vid] && push!(lits, vid + nvars)
        end
        node_lits[node_idx] = lits
    end

    # Compute transitive reachable literals per node via reverse topological sweep.
    reachable_lits = [copy(node_lits[i]) for i in 1:nnodes]
    topo_order = topological_sort_by_dfs(g)
    for node_idx in reverse(topo_order)
        for succ in outneighbors(g, node_idx)
            union!(reachable_lits[node_idx], reachable_lits[succ])
        end
    end

    true_lits = BitSet()

    # If there is a path from ¬v to v, set v to true.
    for vid in 1:nvars
        pos = vid
        neg = vid + nvars
        neg_node = node_for_pos(neg)
        if in(pos, reachable_lits[neg_node])
            pos_node = node_for_pos(pos)
            union!(true_lits, reachable_lits[pos_node])
        end
    end

    # Propagate truth along implications; negations of true literals become false.
    seed_nodes = BitSet(node_for_pos(lit) for lit in true_lits)
    for node_idx in seed_nodes
        union!(true_lits, reachable_lits[node_idx])
    end

    false_lits = BitSet(neg_pos(lit) for lit in true_lits)
    if !isempty(intersect(true_lits, false_lits))
        for lit in true_lits
            lit_is_true[lit] = true
        end
        for lit in false_lits
            lit_is_false[lit] = true
        end
        return false, lit_is_true, lit_is_false
    end

    node_is_true = falses(nnodes)
    node_is_false = falses(nnodes)
    for lit in true_lits
        node_is_true[node_for_pos(lit)] = true
        lit_is_true[lit] = true
    end
    for lit in false_lits
        node_is_false[node_for_pos(lit)] = true
        lit_is_false[lit] = true
    end

    # Infeasible if an implication edge points from true to false.
    for u in 1:nnodes
        node_is_true[u] || continue
        for v in outneighbors(g, u)
            node_is_false[v] && return false, lit_is_true, lit_is_false
        end
    end

    return true, lit_is_true, lit_is_false
end

"""
    get_scc_repr_pos!(prop::ParityPropagator, pos::Int) -> Int

Find and return the union-find representative of literal position `pos`.

Performs path compression on `prop.parent_pos`.
"""
function get_scc_repr_pos!(prop::ParityPropagator, pos::Int)
    # find root
    root = pos
    while prop.parent_pos[root] != root
        root = prop.parent_pos[root]
    end

    # path compression
    while pos != root
        p = prop.parent_pos[pos]
        prop.parent_pos[pos] = root
        pos = p
    end

    return root
end

"""
    union_scc!(prop::ParityPropagator, pos1::Int, pos2::Int) -> Int

Union the SCC representative sets containing `pos1` and `pos2`.

The smaller representative index is kept as root.

# Returns
- The resulting representative literal position.
"""
function union_scc!(prop::ParityPropagator, pos1::Int, pos2::Int)
    root1 = get_scc_repr_pos!(prop, pos1)
    root2 = get_scc_repr_pos!(prop, pos2)
    if root2 < root1
        root1, root2 = root2, root1
    end

    prop.parent_pos[root2] = root1

    return root1
end
