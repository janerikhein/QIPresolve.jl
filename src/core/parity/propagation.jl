using DataStructures: Deque
using Graphs


const VarId = Int

@enum LitLabel TRUE FALSE UNDEF

struct VarLit
    vid::VarId
    neg::Bool
end

"""
    negated(lit::VarLit) -> VarLit

Return the complementary literal of `lit` by flipping its negation flag while
keeping the same variable id.
"""
negated(lit::VarLit) = VarLit(lit.vid, !lit.neg)

mutable struct PropagationManager
    nvars::Int
    nreps::Int

    # Literal index mappings for UnionFind structure
    lit_to_pos::Dict{VarLit, Int}
    pos_to_lit::Vector{VarLit}

    # UnionFind tree structure (parent_pos[i] == i, iff i is representant) 
    parent_pos::Vector{Int}

    # Index mapping for implication graph storing only representant literals
    rep_pos_to_scc_pos::Dict{Int, Int}
    scc_pos_to_rep_pos::Vector{Int}

    # Labels and deques for querying pending fixations and substitutions 
    lit_labels::Vector{LitLabel}
    pending_fixings::Deque{Tuple{VarId, Bool}}
    pending_substitutions::Deque{VarId}
    seen_fixings::Dict{Tuple{VarId, Bool}, Bool}
    seen_substitutions::Dict{VarId, Bool}

    # Implication graph structure over SCC representant positions
    graph::DiGraph 
end

"""
    PropagationManager(var_ids::Vector{VarId}) -> PropagationManager

Construct a propagation manager over the given variable ids.

Each variable is initialized with two literal positions, one for the positive
literal and one for the negated literal. Initially every literal is its own
representative and SCC node, all labels are `UNDEF`, both pending-event queues
are empty, and the implication graph has no edges.
"""
function PropagationManager(var_ids::Vector{VarId})
    nvars = length(var_ids)
    nreps = 2 * nvars

    lit_to_pos = Dict{VarLit, Int}()
    sizehint!(lit_to_pos, nreps)
    pos_to_lit = Vector{VarLit}(undef, nreps)

    @inbounds for (i, vid) in enumerate(var_ids)
        pos = 2 * i - 1
        pos_lit = VarLit(vid, false)
        neg_lit = VarLit(vid, true)
        @assert !haskey(lit_to_pos, pos_lit)
        lit_to_pos[pos_lit] = pos
        lit_to_pos[neg_lit] = pos + 1
        pos_to_lit[pos] = pos_lit
        pos_to_lit[pos + 1] = neg_lit
    end

    parent_pos = collect(1:nreps)
    rep_pos_to_scc_pos = Dict(i => i for i in 1:nreps)
    scc_pos_to_rep_pos = collect(1:nreps)
    lit_labels = fill(UNDEF, nreps)
    pending_fixings = Deque{Tuple{VarId, Bool}}()
    pending_substitutions = Deque{VarId}()
    seen_fixings = Dict{Tuple{VarId, Bool}, Bool}()
    seen_substitutions = Dict{VarId, Bool}()
    graph = DiGraph(nreps)

    return PropagationManager(
        nvars,
        nreps,
        lit_to_pos,
        pos_to_lit,
        parent_pos,
        rep_pos_to_scc_pos,
        scc_pos_to_rep_pos,
        lit_labels,
        pending_fixings,
        pending_substitutions,
        seen_fixings,
        seen_substitutions,
        graph,
    )
end


"""
    repr!(manager::PropagationManager, pos::Int) -> Int

Return the current union-find representative position for `pos`.

This performs path compression on `manager.parent_pos` and therefore mutates
the manager even though the represented partition stays unchanged.
"""
function repr!(manager::PropagationManager, pos::Int)
    root = pos
    while manager.parent_pos[root] != root
        root = manager.parent_pos[root]
    end

    # path compression
    while pos != root
        p = manager.parent_pos[pos]
        manager.parent_pos[pos] = root
        pos = p
    end

    return root
end

"""
    enqueue_fixing!(manager::PropagationManager, vid::VarId, val::Bool) -> Bool

Push `(vid, val)` onto the pending fixings FIFO queue if that exact fixing has
not been queued before.

Returns `true` if a new queue entry was added and `false` if the fixing had
already been recorded earlier.
"""
function enqueue_fixing!(manager::PropagationManager, vid::VarId, val::Bool)
    key = (vid, val)
    get(manager.seen_fixings, key, false) && return false
    manager.seen_fixings[key] = true
    push!(manager.pending_fixings, key)
    return true
end

"""
    enqueue_substitution!(manager::PropagationManager, vid::VarId) -> Bool

Push `vid` onto the pending substitutions FIFO queue if that variable has not
been queued as substitutable before.

Returns `true` if a new queue entry was added and `false` otherwise.
"""
function enqueue_substitution!(manager::PropagationManager, vid::VarId)
    get(manager.seen_substitutions, vid, false) && return false
    manager.seen_substitutions[vid] = true
    push!(manager.pending_substitutions, vid)
    return true
end

"""
    maybe_enqueue_fixing!(manager::PropagationManager, scc_pos::Int)

Inspect the variable represented by SCC position `scc_pos` and enqueue a
variable fixing if the current positive/negative literal labels determine a
truth value.

This is an internal helper used after label updates. If the positive literal is
`TRUE` and the negative literal is `FALSE`, the variable is fixed to `true`. If
the labels are reversed, the variable is fixed to `false`.
"""
function maybe_enqueue_fixing!(manager::PropagationManager, scc_pos::Int)
    lit = manager.pos_to_lit[manager.scc_pos_to_rep_pos[scc_pos]]
    pos_rep = repr!(manager, manager.lit_to_pos[VarLit(lit.vid, false)])
    neg_rep = repr!(manager, manager.lit_to_pos[VarLit(lit.vid, true)])
    pos_scc = manager.rep_pos_to_scc_pos[pos_rep]
    neg_scc = manager.rep_pos_to_scc_pos[neg_rep]
    pos_label = manager.lit_labels[pos_scc]
    neg_label = manager.lit_labels[neg_scc]

    if pos_label == TRUE && neg_label == FALSE
        enqueue_fixing!(manager, lit.vid, true)
    elseif pos_label == FALSE && neg_label == TRUE
        enqueue_fixing!(manager, lit.vid, false)
    end

    return nothing
end

"""
    set_lit_label!(manager::PropagationManager, scc_pos::Int, label::LitLabel) -> Bool

Assign `label` to SCC position `scc_pos` if it is not already set to that
value, then check whether the affected variable has become fixed.

Returns `true` if the label changed and `false` if the existing label already
matched `label`.
"""
function set_lit_label!(manager::PropagationManager, scc_pos::Int, label::LitLabel)
    cur_label = manager.lit_labels[scc_pos]
    cur_label == label && return false
    manager.lit_labels[scc_pos] = label
    maybe_enqueue_fixing!(manager, scc_pos)
    return true
end

"""
    union_reps!(manager::PropagationManager, rep_positions::Vector{Int}) -> Int

Merge the representative positions in `rep_positions` into a single
union-find class and return the chosen representative position.

The minimum position is used as the new representative. Every other merged
representative stops being a representative; if that position corresponds to a
positive literal, its variable is queued for substitution.
"""
function union_reps!(manager::PropagationManager, rep_positions::Vector{Int})
    new_rep = min(rep_positions...)

    for pos in rep_positions
        @assert manager.parent_pos[pos] == pos
        pos == new_rep && continue
        manager.parent_pos[pos] = new_rep
        lit = manager.pos_to_lit[pos]
        !lit.neg && enqueue_substitution!(manager, lit.vid)
    end

    return new_rep
end


"""
    add_edge!(manager::PropagationManager, lit1::VarLit, lit2::VarLit) -> Bool

Add the implication edge `lit1 -> lit2` to the current SCC graph.

The edge is added between the current representative SCC nodes of the two
literals. If both literals already belong to the same representative, no edge
is added. Returns the result of `Graphs.add_edge!`.
"""
function add_edge!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)
    rep1 = repr!(manager, manager.lit_to_pos[lit1])
    rep2 = repr!(manager, manager.lit_to_pos[lit2])
    added = false
    
    if rep1 != rep2
        scc_pos1 = manager.rep_pos_to_scc_pos[rep1]
        scc_pos2 = manager.rep_pos_to_scc_pos[rep2]
        added = Graphs.add_edge!(manager.graph, scc_pos1, scc_pos2)
    end

    return added
end


"""
    add_implication!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)

Insert the implication `lit1 => lit2` into the propagation graph.

To maintain logical equivalence in the implication graph, this also inserts the
contrapositive implication `¬lit2 => ¬lit1`.
"""
function add_implication!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)
    add_edge!(manager, lit1, lit2)
    add_edge!(manager, negated(lit2), negated(lit1))
end


"""
    add_equivalence!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)

Insert a bidirectional equivalence between `lit1` and `lit2`.

This is implemented as the two implications `lit1 => lit2` and `lit2 => lit1`
including their contrapositives.
"""
function add_equivalence!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)
    add_implication!(manager, lit1, lit2)
    add_implication!(manager, lit2, lit1)
end


"""
    fix_var!(manager::PropagationManager, vid::VarId, val::Bool)

Force variable `vid` to take value `val`.

This adds the implication from the opposite literal to the chosen literal,
directly labels the variable's positive and negative literal SCCs accordingly,
and enqueues the fixing in the pending fixings queue if it has not been seen
before.
"""
function fix_var!(manager::PropagationManager, vid::VarId, val::Bool)
    lit = VarLit(vid, !val)
    add_edge!(manager, negated(lit), lit)
    pos_rep = repr!(manager, manager.lit_to_pos[VarLit(vid, false)])
    neg_rep = repr!(manager, manager.lit_to_pos[VarLit(vid, true)])
    pos_scc = manager.rep_pos_to_scc_pos[pos_rep]
    neg_scc = manager.rep_pos_to_scc_pos[neg_rep]
    if val
        set_lit_label!(manager, pos_scc, TRUE)
        set_lit_label!(manager, neg_scc, FALSE)
    else
        set_lit_label!(manager, pos_scc, FALSE)
        set_lit_label!(manager, neg_scc, TRUE)
    end
    enqueue_fixing!(manager, vid, val)
end


"""
    update_sccs!(manager::PropagationManager) -> PropagationManager

Recompute strongly connected components of the current implication graph and
rebuild the condensed SCC graph.

This updates union-find representatives, SCC index mappings, the condensed graph
itself, and the SCC-level label vector. Labels from merged SCC nodes are carried
over to the new SCC representatives, and substitution events are queued through
`union_reps!`.
"""
function update_sccs!(manager::PropagationManager)
    old_graph = manager.graph
    old_scc_pos_to_rep_pos = manager.scc_pos_to_rep_pos
    old_lit_labels = manager.lit_labels

    comps = strongly_connected_components(old_graph)
    ncomps = length(comps)

    new_rep_to_idx = Dict{Int, Int}()
    sizehint!(new_rep_to_idx, ncomps)
    new_idx_to_rep = Vector{Int}(undef, ncomps)
    new_lit_labels = fill(UNDEF, ncomps)

    @inbounds for (i, comp_scc_pos) in enumerate(comps)
        comp_rep_pos = old_scc_pos_to_rep_pos[comp_scc_pos]
        new_rep = union_reps!(manager, comp_rep_pos)
        new_rep_to_idx[new_rep] = i
        new_idx_to_rep[i] = new_rep

        label = UNDEF
        for v in comp_scc_pos
            cur_label = old_lit_labels[v]
            cur_label == UNDEF && continue
            label = cur_label
            break
        end
        new_lit_labels[i] = label
    end

    new_graph = DiGraph(ncomps)
    @inbounds for edge in edges(old_graph)
        src_rep = repr!(manager, old_scc_pos_to_rep_pos[src(edge)])
        dst_rep = repr!(manager, old_scc_pos_to_rep_pos[dst(edge)])
        if src_rep != dst_rep
            src_idx = new_rep_to_idx[src_rep]
            dst_idx = new_rep_to_idx[dst_rep]
            Graphs.add_edge!(new_graph, src_idx, dst_idx)
        end
    end

    manager.nreps = ncomps
    manager.rep_pos_to_scc_pos = new_rep_to_idx
    manager.scc_pos_to_rep_pos = new_idx_to_rep
    manager.lit_labels = new_lit_labels
    manager.graph = new_graph
    return manager
end

"""
    propagate_labels!(manager::PropagationManager, top_order) -> PropagationManager

Propagate existing SCC labels through the DAG defined by `top_order`.

`TRUE` labels are propagated forward along outgoing implications, and `FALSE`
labels are propagated backward along incoming implications. Newly assigned SCC
labels are routed through `set_lit_label!` so any implied variable fixings are
queued immediately.
"""
function propagate_labels!(manager::PropagationManager, top_order)
    graph = manager.graph
    lit_labels = manager.lit_labels

    # Forward: TRUE propagates along outgoing implications.
    @inbounds for v in top_order
        lit_labels[v] == TRUE || continue
        for w in outneighbors(graph, v)
            lit_labels[w] == UNDEF && set_lit_label!(manager, w, TRUE)
        end
    end

    # Backward: FALSE propagates along incoming implications.
    @inbounds for i in length(top_order):-1:1
        w = top_order[i]
        lit_labels[w] == FALSE || continue
        for v in inneighbors(graph, w)
            lit_labels[v] == UNDEF && set_lit_label!(manager, v, FALSE)
        end
    end

    return manager
end


"""
    label_reachables!(manager::PropagationManager, top_order) -> PropagationManager

Detect forced literals via reachability in the condensed implication graph.

For each SCC node `v`, this computes whether the SCC of the opposite literal can
reach `v`. If so, the opposite SCC is labeled `FALSE` and `v` is labeled
`TRUE`, again using `set_lit_label!` so any resulting variable fixings are
queued.
"""
function label_reachables!(manager::PropagationManager, top_order)
    n = nv(manager.graph)
    reachables = [falses(n) for _ in 1:n]
    neg_nodes = Vector{Int}(undef, n)

    @inbounds for v in 1:n
        lit = manager.pos_to_lit[manager.scc_pos_to_rep_pos[v]]
        neg_rep = repr!(manager, manager.lit_to_pos[negated(lit)])
        neg_nodes[v] = manager.rep_pos_to_scc_pos[neg_rep]
    end

    @inbounds for i in length(top_order):-1:1
        v = top_order[i]
        reachable_v = reachables[v]
        for w in outneighbors(manager.graph, v)
            reachable_v[w] = true
            reachable_v .|= reachables[w]
        end
    end

    @inbounds for v in top_order
        neg_v = neg_nodes[v]
        if reachables[neg_v][v]
            set_lit_label!(manager, neg_v, FALSE)
            set_lit_label!(manager, v, TRUE)
        end
    end

    return manager
end

"""
    pop_fixing!(manager::PropagationManager) -> Union{Nothing, Tuple{VarId, Bool}}

Pop and return the oldest pending variable fixing from the FIFO queue.

Returns `nothing` if no pending fixing is available.
"""
function pop_fixing!(manager::PropagationManager)
    isempty(manager.pending_fixings) && return nothing
    return popfirst!(manager.pending_fixings)
end

"""
    pop_substitution!(manager::PropagationManager) -> Union{Nothing, Tuple{VarId, VarId, Bool}}

Pop and return the oldest pending variable substitution from the FIFO queue.

The returned tuple is `(var_id, rep_var_id, negated)`, where `rep_var_id` is
the current representative variable reached from the positive literal of
`var_id`, and `negated` indicates whether that representative literal is
negated. Returns `nothing` if the queue is empty.
"""
function pop_substitution!(manager::PropagationManager)
    isempty(manager.pending_substitutions) && return nothing
    var_id = popfirst!(manager.pending_substitutions)
    rep_pos = repr!(manager, manager.lit_to_pos[VarLit(var_id, false)])
    rep_lit = manager.pos_to_lit[rep_pos]
    return (var_id, rep_lit.vid, rep_lit.neg)
end

"""
    update!(manager::PropagationManager) -> PropagationManager

Refresh the propagation state after graph changes.

This recomputes SCCs, builds a topological order of the condensed graph, labels
forced literals via reachability, and then propagates those labels through the
implication DAG.
"""
function update!(manager::PropagationManager)
    update_sccs!(manager)
    top_order = topological_sort(manager.graph)
    label_reachables!(manager, top_order)
    propagate_labels!(manager, top_order)
    return manager
end


function register_implications!(manager::PropagationManager, con::XorConstraint, idx_to_vid)
    con.meta.requires_prop || return nothing
    con.meta.requires_update && update!(con)

    # Case 1: Pure XOR
    if con.meta.nnz_conj == 0 && con.meta.nnz_par > 0

        # Case 1.1: singleton parity
        if con.meta.nnz_par == 1
            vid = idx_to_vid[findfirst(con.par)]
            fix_var!(manager, vid, con.rhs)
        end

        # Case 1.2: two parity terms
        if con.meta.nnz_par == 2
            vid1_idx = findfirst(con.par)
            vid1 = idx_to_vid[vid1_idx]
            vid2 = idx_to_vid[findnext(con.par, vid1_idx + 1)]
            lit1 = VarLit(vid1, false)
            lit2 = VarLit(vid2, con.rhs)
            add_equivalence!(manager, lit1, lit2)
        end
    end

    # Case 2: Pure conjunction
    if con.meta.nnz_conj > 0 && con.meta.nnz_par == 0

        # Case 2.1: singleton conjunction term
        if con.meta.nnz_conj == 1
            (vid1_idx, vid2_idx) = findfirst(con.conj).I
            vid1, vid2 = idx_to_vid[vid1_idx], idx_to_vid[vid2_idx]
            if con.rhs 
                fix_var!(manager, vid1, true)
                fix_var!(manager, vid2, true)
            else
                lit1 = VarLit(vid1, false)
                lit2 = VarLit(vid2, true)
                add_implication!(manager, lit1, lit2)
            end
        end
    end

    con.meta.requires_prop = false

end
