using DataStructures: Deque
using Graphs

const VarId = Int

"""
    LitLabel

Represent the propagated truth status of one SCC literal node.
"""
@enum LitLabel TRUE FALSE UNDEF

"""
    VarLit

Represent one Boolean literal associated with a model variable.

# Fields
- `vid`: Variable id.
- `neg`: Whether the literal is negated.
"""
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

"""
    PropagationManager

Store the union-find and implication-graph state used for parity propagation.

This structure tracks literal representatives, SCC labels, pending fixings and
substitutions, and the condensed implication graph across propagation passes.
"""
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
    ParityPropagator

Alias `PropagationManager` for parity-specific call sites.
"""
const ParityPropagator = PropagationManager

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
    ensure_literals!(manager::PropagationManager, var_ids::Vector{VarId}) -> PropagationManager

Ensure that every variable id in `var_ids` has positive and negated literals
registered in `manager`.

Missing variables are appended as isolated singleton SCCs with `UNDEF` labels.
Existing literals are left unchanged.
"""
function ensure_literals!(manager::PropagationManager, var_ids::Vector{VarId})
    for vid in var_ids
        haskey(manager.lit_to_pos, VarLit(vid, false)) && continue

        pos_pos = length(manager.pos_to_lit) + 1
        neg_pos = pos_pos + 1

        push!(manager.pos_to_lit, VarLit(vid, false))
        push!(manager.pos_to_lit, VarLit(vid, true))
        manager.lit_to_pos[VarLit(vid, false)] = pos_pos
        manager.lit_to_pos[VarLit(vid, true)] = neg_pos

        push!(manager.parent_pos, pos_pos)
        push!(manager.parent_pos, neg_pos)

        pos_scc = manager.nreps + 1
        neg_scc = manager.nreps + 2
        manager.rep_pos_to_scc_pos[pos_pos] = pos_scc
        manager.rep_pos_to_scc_pos[neg_pos] = neg_scc
        push!(manager.scc_pos_to_rep_pos, pos_pos)
        push!(manager.scc_pos_to_rep_pos, neg_pos)
        push!(manager.lit_labels, UNDEF)
        push!(manager.lit_labels, UNDEF)
        Graphs.add_vertex!(manager.graph)
        Graphs.add_vertex!(manager.graph)

        manager.nvars += 1
        manager.nreps += 2
    end

    return manager
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
function enqueue_fixing!(
    manager::PropagationManager,
    vid::VarId,
    val::Bool,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
    source::Symbol = :unknown,
)
    key = (vid, val)
    get(manager.seen_fixings, key, false) && return false
    manager.seen_fixings[key] = true
    push!(manager.pending_fixings, key)
    _record_parity_fixing!(stats_accumulator, source)
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
function maybe_enqueue_fixing!(
    manager::PropagationManager,
    scc_pos::Int,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
    source::Symbol = :unknown,
)
    lit = manager.pos_to_lit[manager.scc_pos_to_rep_pos[scc_pos]]
    pos_rep = repr!(manager, manager.lit_to_pos[VarLit(lit.vid, false)])
    neg_rep = repr!(manager, manager.lit_to_pos[VarLit(lit.vid, true)])
    pos_scc = manager.rep_pos_to_scc_pos[pos_rep]
    neg_scc = manager.rep_pos_to_scc_pos[neg_rep]
    pos_label = manager.lit_labels[pos_scc]
    neg_label = manager.lit_labels[neg_scc]

    if pos_label == TRUE && neg_label == FALSE
        enqueue_fixing!(manager, lit.vid, true, stats_accumulator, source)
    elseif pos_label == FALSE && neg_label == TRUE
        enqueue_fixing!(manager, lit.vid, false, stats_accumulator, source)
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
function set_lit_label!(
    manager::PropagationManager,
    scc_pos::Int,
    label::LitLabel,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
    source::Symbol = :unknown,
)
    cur_label = manager.lit_labels[scc_pos]
    cur_label == label && return false
    manager.lit_labels[scc_pos] = label
    maybe_enqueue_fixing!(manager, scc_pos, stats_accumulator, source)
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
function add_edge!(
    manager::PropagationManager,
    lit1::VarLit,
    lit2::VarLit,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
)
    rep1 = repr!(manager, manager.lit_to_pos[lit1])
    rep2 = repr!(manager, manager.lit_to_pos[lit2])
    added = false
    
    if rep1 != rep2
        scc_pos1 = manager.rep_pos_to_scc_pos[rep1]
        scc_pos2 = manager.rep_pos_to_scc_pos[rep2]
        added = Graphs.add_edge!(manager.graph, scc_pos1, scc_pos2)
    end

    _record_implication_edge!(stats_accumulator, added)
    return added
end


"""
    add_implication!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)

Insert the implication `lit1 => lit2` into the propagation graph.

To maintain logical equivalence in the implication graph, this also inserts the
contrapositive implication `¬lit2 => ¬lit1`.
"""
function add_implication!(
    manager::PropagationManager,
    lit1::VarLit,
    lit2::VarLit,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
)
    add_edge!(manager, lit1, lit2, stats_accumulator)
    add_edge!(manager, negated(lit2), negated(lit1), stats_accumulator)
end


"""
    add_equivalence!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)

Insert a bidirectional equivalence between `lit1` and `lit2`.

This is implemented as the two implications `lit1 => lit2` and `lit2 => lit1`
including their contrapositives.
"""
function add_equivalence!(
    manager::PropagationManager,
    lit1::VarLit,
    lit2::VarLit,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
)
    add_implication!(manager, lit1, lit2, stats_accumulator)
    add_implication!(manager, lit2, lit1, stats_accumulator)
end


"""
    fix_var!(manager::PropagationManager, vid::VarId, val::Bool)

Force variable `vid` to take value `val`.

This adds the implication from the opposite literal to the chosen literal,
directly labels the variable's positive and negative literal SCCs accordingly,
and enqueues the fixing in the pending fixings queue if it has not been seen
before.
"""
function fix_var!(
    manager::PropagationManager,
    vid::VarId,
    val::Bool,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
    source::Symbol = :elimination,
)
    lit = VarLit(vid, !val)
    add_edge!(manager, negated(lit), lit, stats_accumulator)
    pos_rep = repr!(manager, manager.lit_to_pos[VarLit(vid, false)])
    neg_rep = repr!(manager, manager.lit_to_pos[VarLit(vid, true)])
    pos_scc = manager.rep_pos_to_scc_pos[pos_rep]
    neg_scc = manager.rep_pos_to_scc_pos[neg_rep]
    if val
        set_lit_label!(manager, pos_scc, TRUE, stats_accumulator, source)
        set_lit_label!(manager, neg_scc, FALSE, stats_accumulator, source)
    else
        set_lit_label!(manager, pos_scc, FALSE, stats_accumulator, source)
        set_lit_label!(manager, neg_scc, TRUE, stats_accumulator, source)
    end
    enqueue_fixing!(manager, vid, val, stats_accumulator, source)
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
function propagate_labels!(
    manager::PropagationManager,
    top_order,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
)
    graph = manager.graph
    lit_labels = manager.lit_labels

    # Forward: TRUE propagates along outgoing implications.
    @inbounds for v in top_order
        lit_labels[v] == TRUE || continue
        for w in outneighbors(graph, v)
            lit_labels[w] == UNDEF && set_lit_label!(manager, w, TRUE, stats_accumulator, :propagation)
        end
    end

    # Backward: FALSE propagates along incoming implications.
    @inbounds for i in length(top_order):-1:1
        w = top_order[i]
        lit_labels[w] == FALSE || continue
        for v in inneighbors(graph, w)
            lit_labels[v] == UNDEF && set_lit_label!(manager, v, FALSE, stats_accumulator, :propagation)
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
function label_reachables!(
    manager::PropagationManager,
    top_order,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
)
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
            set_lit_label!(manager, neg_v, FALSE, stats_accumulator, :path)
            set_lit_label!(manager, v, TRUE, stats_accumulator, :path)
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
    fixed_value(manager::PropagationManager, vid::VarId) -> Union{Nothing, Bool}

Return the fixed Boolean value implied for `vid`, if any.

This inspects the current labels of the positive and negated literal SCCs for
`vid`. Returns `true` when `vid` is fixed true, `false` when fixed false, and
`nothing` when no fixing is currently implied or `vid` is unknown to the
manager.
"""
function fixed_value(manager::PropagationManager, vid::VarId)
    pos_pos = get(manager.lit_to_pos, VarLit(vid, false), 0)
    pos_pos == 0 && return nothing

    neg_pos = manager.lit_to_pos[VarLit(vid, true)]
    pos_rep = repr!(manager, pos_pos)
    neg_rep = repr!(manager, neg_pos)
    pos_scc = manager.rep_pos_to_scc_pos[pos_rep]
    neg_scc = manager.rep_pos_to_scc_pos[neg_rep]
    pos_label = manager.lit_labels[pos_scc]
    neg_label = manager.lit_labels[neg_scc]

    if pos_label == TRUE && neg_label == FALSE
        return true
    elseif pos_label == FALSE && neg_label == TRUE
        return false
    end

    return nothing
end

"""
    _component_size(manager, rep_pos)

Count how many literal positions currently belong to representative `rep_pos`.
"""
function _component_size(manager::PropagationManager, rep_pos::Int)
    return count(pos -> repr!(manager, pos) == rep_pos, eachindex(manager.pos_to_lit))
end

"""
    _rep_component_positions(manager, rep_pos)

Collect the literal positions currently represented by `rep_pos`.
"""
function _rep_component_positions(manager::PropagationManager, rep_pos::Int)
    return [pos for pos in eachindex(manager.pos_to_lit) if repr!(manager, pos) == rep_pos]
end

"""
    substitute_scc_by_new_var!(
        manager::PropagationManager,
        scc_vid::VarId,
        new_vid::VarId,
    ) -> PropagationManager

Replace the current positive and negative SCCs of `scc_vid` by fresh literals
for `new_vid`.

The current SCC ids of `scc_vid` are reused for `new_vid`. Every old literal in
the substituted positive and negative components is split back into its own
singleton SCC, and those fresh SCC nodes are left isolated in the graph.

The affected SCCs must both be unlabeled, have size at least two, and the
manager must not already contain `new_vid`.
"""
function substitute_scc_by_new_var!(
    manager::PropagationManager,
    scc_vid::VarId,
    new_vid::VarId,
)
    @assert !haskey(manager.lit_to_pos, VarLit(new_vid, false))
    @assert !haskey(manager.lit_to_pos, VarLit(new_vid, true))
    @assert isempty(manager.pending_substitutions)

    pos_pos = manager.lit_to_pos[VarLit(scc_vid, false)]
    neg_pos = manager.lit_to_pos[VarLit(scc_vid, true)]
    pos_rep = repr!(manager, pos_pos)
    neg_rep = repr!(manager, neg_pos)
    @assert pos_rep != neg_rep

    pos_scc = manager.rep_pos_to_scc_pos[pos_rep]
    neg_scc = manager.rep_pos_to_scc_pos[neg_rep]
    @assert pos_scc != neg_scc
    @assert manager.lit_labels[pos_scc] == UNDEF
    @assert manager.lit_labels[neg_scc] == UNDEF
    pos_component = [pos for pos in eachindex(manager.pos_to_lit) if repr!(manager, pos) == pos_rep]
    neg_component = [pos for pos in eachindex(manager.pos_to_lit) if repr!(manager, pos) == neg_rep]
    @assert length(pos_component) >= 2
    @assert length(neg_component) >= 2

    new_pos_rep = length(manager.pos_to_lit) + 1
    new_neg_rep = new_pos_rep + 1
    old_nreps = manager.nreps
    split_positions = vcat(pos_component, neg_component)

    push!(manager.pos_to_lit, VarLit(new_vid, false))
    push!(manager.pos_to_lit, VarLit(new_vid, true))
    manager.lit_to_pos[VarLit(new_vid, false)] = new_pos_rep
    manager.lit_to_pos[VarLit(new_vid, true)] = new_neg_rep
    push!(manager.parent_pos, new_pos_rep)
    push!(manager.parent_pos, new_neg_rep)

    for pos in split_positions
        manager.parent_pos[pos] = pos
    end

    manager.nvars += 1
    manager.rep_pos_to_scc_pos[new_pos_rep] = pos_scc
    manager.rep_pos_to_scc_pos[new_neg_rep] = neg_scc
    manager.scc_pos_to_rep_pos[pos_scc] = new_pos_rep
    manager.scc_pos_to_rep_pos[neg_scc] = new_neg_rep

    next_scc = old_nreps + 1
    for pos in split_positions
        manager.rep_pos_to_scc_pos[pos] = next_scc
        push!(manager.scc_pos_to_rep_pos, pos)
        push!(manager.lit_labels, UNDEF)
        Graphs.add_vertex!(manager.graph)
        next_scc += 1
    end

    manager.nreps = next_scc - 1

    return manager
end

"""
    finalize_phase!(manager::PropagationManager) -> PropagationManager

Prepare `manager` for the next parity presolve phase while preserving all
unlabeled SCC structure.

Every labeled SCC is stripped of incident graph edges, all literals in the
labeled component become singleton SCCs, and all affected SCC labels are reset
to `UNDEF`.
"""
function finalize_phase!(manager::PropagationManager)
    @assert isempty(manager.pending_fixings)
    @assert isempty(manager.pending_substitutions)

    labeled_sccs = [scc for scc in eachindex(manager.lit_labels) if manager.lit_labels[scc] != UNDEF]
    isempty(labeled_sccs) || begin
        labeled_scc_set = Set(labeled_sccs)
        old_graph = manager.graph
        new_graph = DiGraph(nv(old_graph))

        for edge in edges(old_graph)
            src(edge) in labeled_scc_set && continue
            dst(edge) in labeled_scc_set && continue
            Graphs.add_edge!(new_graph, src(edge), dst(edge))
        end
        manager.graph = new_graph

        next_scc = manager.nreps + 1
        for scc in labeled_sccs
            rep_pos = manager.scc_pos_to_rep_pos[scc]
            positions = _rep_component_positions(manager, rep_pos)

            for pos in positions
                manager.parent_pos[pos] = pos
            end

            manager.rep_pos_to_scc_pos[rep_pos] = scc
            manager.scc_pos_to_rep_pos[scc] = rep_pos
            manager.lit_labels[scc] = UNDEF

            for pos in positions
                pos == rep_pos && continue
                manager.rep_pos_to_scc_pos[pos] = next_scc
                push!(manager.scc_pos_to_rep_pos, pos)
                push!(manager.lit_labels, UNDEF)
                Graphs.add_vertex!(manager.graph)
                next_scc += 1
            end
        end

        manager.nreps = next_scc - 1
    end

    empty!(manager.seen_fixings)
    empty!(manager.seen_substitutions)
    return manager
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
function update!(
    manager::PropagationManager,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
)
    update_sccs!(manager)
    top_order = topological_sort(manager.graph)
    label_reachables!(manager, top_order, stats_accumulator)
    propagate_labels!(manager, top_order, stats_accumulator)
    return manager
end

@inline _lit_for_value(vid::VarId, val::Bool) = VarLit(vid, !val)

@inline _lit_sort_key(lit::VarLit) = (lit.vid, lit.neg)

function _canonical_implication_key(antecedent::VarLit, consequent::VarLit)
    key = (_lit_sort_key(antecedent), _lit_sort_key(consequent))
    contra_key = (_lit_sort_key(negated(consequent)), _lit_sort_key(negated(antecedent)))
    return isless(contra_key, key) ? contra_key : key
end

function _active_conjunction_pairs(conj::Union{BitMatrix, Nothing})
    conj === nothing && return Tuple{Int, Int}[]

    pairs = Tuple{Int, Int}[]
    @inbounds for row in 1:(size(conj, 1) - 1)
        for col in (row + 1):size(conj, 1)
            conj[row, col] && push!(pairs, (row, col))
        end
    end
    return pairs
end

function _two_term_satisfying_assignments(
        con::XorConstraint,
        parity_indices::Vector{Int},
        conjunction_pairs::Vector{Tuple{Int, Int}},
        involved_indices::Vector{Int},
    )
    index_to_position = Dict(idx => pos for (pos, idx) in enumerate(involved_indices))
    assignments = BitVector[]
    n = length(involved_indices)

    for mask in 0:(2^n - 1)
        assignment = falses(n)
        @inbounds for pos in 1:n
            assignment[pos] = ((mask >> (pos - 1)) & 1) == 1
        end

        value = false
        for idx in parity_indices
            value ⊻= assignment[index_to_position[idx]]
        end
        for (first_idx, second_idx) in conjunction_pairs
            value ⊻= assignment[index_to_position[first_idx]] &
                assignment[index_to_position[second_idx]]
        end

        value == con.rhs && push!(assignments, assignment)
    end

    return assignments
end

function _add_deduplicated_implication!(
        manager::PropagationManager,
        antecedent::VarLit,
        consequent::VarLit,
        seen_implications::Set{Tuple{Tuple{VarId, Bool}, Tuple{VarId, Bool}}},
        stats_accumulator::Union{Nothing, _ParityStatsAccumulator},
    )
    antecedent == consequent && return manager

    key = _canonical_implication_key(antecedent, consequent)
    key in seen_implications && return manager
    push!(seen_implications, key)
    add_implication!(manager, antecedent, consequent, stats_accumulator)
    return manager
end

function _register_two_term_xor_and_implications!(
        manager::PropagationManager,
        con::XorConstraint,
        idx_to_vid,
        stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
    )
    parity_indices = findall(con.par)
    conjunction_pairs = _active_conjunction_pairs(con.conj)
    length(parity_indices) + length(conjunction_pairs) == 2 || return false
    isempty(conjunction_pairs) && return false

    involved_indices = Int[]
    append!(involved_indices, parity_indices)
    for (first_idx, second_idx) in conjunction_pairs
        push!(involved_indices, first_idx)
        push!(involved_indices, second_idx)
    end
    sort!(unique!(involved_indices))

    satisfying_assignments = _two_term_satisfying_assignments(
        con,
        parity_indices,
        conjunction_pairs,
        involved_indices,
    )
    isempty(satisfying_assignments) && return true

    seen_implications = Set{Tuple{Tuple{VarId, Bool}, Tuple{VarId, Bool}}}()
    for (antecedent_pos, antecedent_idx) in enumerate(involved_indices)
        antecedent_vid = idx_to_vid[antecedent_idx]
        for antecedent_value in (false, true)
            matching_assignments = [
                assignment for assignment in satisfying_assignments
                if assignment[antecedent_pos] == antecedent_value
            ]

            if isempty(matching_assignments)
                _add_deduplicated_implication!(
                    manager,
                    _lit_for_value(antecedent_vid, antecedent_value),
                    _lit_for_value(antecedent_vid, !antecedent_value),
                    seen_implications,
                    stats_accumulator,
                )
                continue
            end

            for (consequent_pos, consequent_idx) in enumerate(involved_indices)
                consequent_pos == antecedent_pos && continue

                consequent_value = matching_assignments[1][consequent_pos]
                all(
                    assignment -> assignment[consequent_pos] == consequent_value,
                    matching_assignments,
                ) || continue

                _add_deduplicated_implication!(
                    manager,
                    _lit_for_value(antecedent_vid, antecedent_value),
                    _lit_for_value(idx_to_vid[consequent_idx], consequent_value),
                    seen_implications,
                    stats_accumulator,
                )
            end
        end
    end

    return true
end

"""
    _triangle_column_matches(conj, col, neighbor1, neighbor2)

Check whether column `col` has exactly the two expected triangle neighbors.
"""
function _triangle_column_matches(
    conj::BitMatrix,
    col::Int,
    neighbor1::Int,
    neighbor2::Int,
)
    count = 0
    n = size(conj, 1)

    @inbounds for row in 1:n
        conj[row, col] || continue
        (row == neighbor1 || row == neighbor2) || return false
        count += 1
    end

    return count == 2
end

"""
    _triangle_vertices(conj)

Detect whether `conj` contains exactly one triangle support and return its
vertex positions.

# Returns
- A tuple `(i, j, k)` when `conj` contains a triangle on exactly three
  vertices.
- `nothing` otherwise.
"""
function _triangle_vertices(conj::BitMatrix)
    n = size(conj, 1)

    @inbounds for col in 1:n
        neighbor1 = findfirst(@view conj[:, col])
        neighbor1 === nothing && continue

        count = 0
        neighbor2 = 0

        for row in 1:n
            conj[row, col] || continue
            count += 1

            if count == 1
                row == neighbor1 || return nothing
            elseif count == 2
                neighbor2 = row
            else
                return nothing
            end
        end

        count == 2 || return nothing
        _triangle_column_matches(conj, col, neighbor1, neighbor2) || return nothing
        _triangle_column_matches(conj, neighbor1, col, neighbor2) || return nothing
        _triangle_column_matches(conj, neighbor2, col, neighbor1) || return nothing
        return (col, neighbor1, neighbor2)
    end

    return nothing
end


"""
    register_implications!(manager, con, idx_to_vid)

Register the direct logical consequences of `con` in `manager`.

Handle supported small pure-XOR, pure-conjunction, and triangle-shaped rows and
mark `con` as no longer requiring propagation afterward.

# Arguments
- `manager`: Propagation manager mutated in place.
- `con`: Constraint whose implications should be registered.
- `idx_to_vid`: Mapping from row indices to model variable ids.

# Returns
- `nothing`.
"""
function register_implications!(
    manager::PropagationManager,
    con::XorConstraint,
    idx_to_vid,
    stats_accumulator::Union{Nothing, _ParityStatsAccumulator} = nothing,
)
    con.meta.requires_prop || return nothing
    con.meta.requires_update && update!(con)

    handled_two_term_xor_and = _register_two_term_xor_and_implications!(
        manager,
        con,
        idx_to_vid,
        stats_accumulator,
    )

    # Case 1: Pure XOR
    if !handled_two_term_xor_and && con.meta.nnz_conj == 0 && con.meta.nnz_par > 0

        # Case 1.1: singleton parity
        if con.meta.nnz_par == 1
            vid = idx_to_vid[findfirst(con.par)]
            fix_var!(manager, vid, con.rhs, stats_accumulator, :elimination)
        end

        # Case 1.2: two parity terms
        if con.meta.nnz_par == 2
            vid1_idx = findfirst(con.par)
            vid1 = idx_to_vid[vid1_idx]
            vid2 = idx_to_vid[findnext(con.par, vid1_idx + 1)]
            lit1 = VarLit(vid1, false)
            lit2 = VarLit(vid2, con.rhs)
            add_equivalence!(manager, lit1, lit2, stats_accumulator)
        end
    end

    # Case 2: Pure conjunction
    if !handled_two_term_xor_and && con.meta.nnz_conj > 0 && con.meta.nnz_par == 0
        conj = con.conj::BitMatrix

        # Case 2.1: singleton conjunction term
        if con.meta.nnz_conj == 1
            vid1_idx = vid2_idx = 0
            @inbounds for i in 1:(size(conj, 1) - 1)
                for j in (i + 1):size(conj, 1)
                    conj[i, j] || continue
                    vid1_idx = i
                    vid2_idx = j
                    break
                end
                vid1_idx != 0 && break
            end

            @assert vid1_idx != 0
            vid1, vid2 = idx_to_vid[vid1_idx], idx_to_vid[vid2_idx]
            if con.rhs 
                fix_var!(manager, vid1, true, stats_accumulator, :elimination)
                fix_var!(manager, vid2, true, stats_accumulator, :elimination)
            else
                lit1 = VarLit(vid1, false)
                lit2 = VarLit(vid2, true)
                add_implication!(manager, lit1, lit2, stats_accumulator)
            end
        end

        # Case 2.2: two conjunction terms
        if con.meta.nnz_conj == 2 && con.rhs
            i = j = k = l = 0
            found = 0

            @inbounds for row in 1:(size(conj, 1) - 1)
                for col in (row + 1):size(conj, 1)
                    conj[row, col] || continue
                    found += 1

                    if found == 1
                        i, j = row, col
                    else
                        @assert found == 2
                        k, l = row, col
                    end
                end
            end

            @assert found == 2
            vids1 = (idx_to_vid[i], idx_to_vid[j])
            vids2 = (idx_to_vid[k], idx_to_vid[l])

            for antecedent_vid in vids1, consequent_vid in vids2
                add_implication!(
                    manager,
                    _lit_for_value(antecedent_vid, false),
                    _lit_for_value(consequent_vid, true),
                    stats_accumulator,
                )
            end

            for antecedent_vid in vids2, consequent_vid in vids1
                add_implication!(
                    manager,
                    _lit_for_value(antecedent_vid, false),
                    _lit_for_value(consequent_vid, true),
                    stats_accumulator,
                )
            end
        end

        # Case 2.3: triangle-shaped conjunction terms
        if con.meta.nnz_conj == 3
            triangle = _triangle_vertices(conj)

            if triangle !== nothing
                antecedent_val = !con.rhs
                consequent_val = con.rhs
                vids = (
                    idx_to_vid[triangle[1]],
                    idx_to_vid[triangle[2]],
                    idx_to_vid[triangle[3]],
                )

                for i in eachindex(vids)
                    antecedent_vid = vids[i]
                    for j in eachindex(vids)
                        i == j && continue
                        consequent_vid = vids[j]
                        add_implication!(
                            manager,
                            _lit_for_value(antecedent_vid, antecedent_val),
                            _lit_for_value(consequent_vid, consequent_val),
                            stats_accumulator,
                        )
                    end
                end
            end
        end
    end

    con.meta.requires_prop = false

end
