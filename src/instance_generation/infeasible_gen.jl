using Graphs
using JuMP
using Random


const BoundingBoxConflict = NamedTuple{
    (:variable, :original_lower, :original_upper, :proposed_lower, :proposed_upper),
    Tuple{String, Float64, Float64, Float64, Float64},
}


"""
Raised when a proposed bounding-box restriction has an empty intersection with
one or more existing variable domains.
"""
struct BoundingBoxCandidateRejected <: Exception
    conflicts::Vector{BoundingBoxConflict}
end


function Base.showerror(io::IO, error::BoundingBoxCandidateRejected)
    print(io, "bounding-box candidate creates empty variable domains: ")
    for (index, conflict) in enumerate(error.conflicts)
        index > 1 && print(io, "; ")
        print(
            io,
            conflict.variable,
            " original=[",
            conflict.original_lower,
            ", ",
            conflict.original_upper,
            "] proposed=[",
            conflict.proposed_lower,
            ", ",
            conflict.proposed_upper,
            "]",
        )
    end
    return nothing
end


function validate_box_margin(box_margin::Integer)::Int
    box_margin isa Bool && throw(ArgumentError(
        "box_margin must be a positive integer, got $box_margin"
    ))
    margin_int = Int(box_margin)
    margin_int >= 1 || throw(ArgumentError(
        "box_margin must be a positive integer, got $box_margin"
    ))
    return margin_int
end


function validate_infeasible_strategy(strategy::Symbol)::Symbol
    if strategy == :bounding_box
        return strategy
    elseif strategy == :vertex_contraction || strategy == :contraction
        return :vertex_contraction
    end
    throw(ArgumentError(
        "strategy must be one of :bounding_box, :vertex_contraction, or :contraction, got $strategy"
    ))
end


function validate_infeasible_base(base::Symbol)::Symbol
    if base == :globally_rigid || base == :laman || base == :random_2_connected
        return base
    elseif base == :two_connected
        return :random_2_connected
    end
    throw(ArgumentError(
        "base must be one of :globally_rigid, :laman, :random_2_connected, or :two_connected, got $base"
    ))
end


function _sample_infeasible_base(
        rng::AbstractRNG, n::Int, base::Symbol; kwargs...
    )::EmbeddedGraph
    normalized_base = validate_infeasible_base(base)
    graph, coords = if normalized_base == :globally_rigid
        _random_globally_rigid_graph(rng, n; kwargs...)
    elseif normalized_base == :laman
        _random_laman_graph(rng, n; kwargs...)
    else
        _random_2_connected_graph(rng, n; kwargs...)
    end
    return to_embedded(graph, coords)
end


function _graph_center_vertex(emb_graph::EmbeddedGraph)::Int
    pairwise_distances = pw_shortest_paths(emb_graph.graph)
    center = graph_center(pairwise_distances)
    isfinite(maximum(@view pairwise_distances[center, :])) ||
        throw(ArgumentError("embedding graph must be connected"))
    return center
end


function _model_frame_coords(
        emb_graph::EmbeddedGraph, center::Int, anchors::AbstractVector{Int}
    )::Vector{IPoint}
    if isempty(anchors)
        center_point = emb_graph.coords[center]
        return [
            IPoint(point.x - center_point.x, point.y - center_point.y)
            for point in emb_graph.coords
        ]
    end
    return copy(emb_graph.coords)
end


function _tightened_bounding_box(
        coords::AbstractVector{IPoint}, box_margin::Int
    )::NTuple{4, Int}
    x_min = minimum(point.x for point in coords)
    x_max = maximum(point.x for point in coords)
    y_min = minimum(point.y for point in coords)
    y_max = maximum(point.y for point in coords)

    x_lower = x_min + box_margin
    x_upper = x_max - box_margin
    y_lower = y_min + box_margin
    y_upper = y_max - box_margin
    return x_lower, x_upper, y_lower, y_upper
end


function _proposed_variable_bounds(var, lower::Int, upper::Int)
    original_lower = lower_bound(var)
    original_upper = upper_bound(var)
    return (
        original_lower = original_lower,
        original_upper = original_upper,
        proposed_lower = max(original_lower, lower),
        proposed_upper = min(original_upper, upper),
    )
end


function _apply_bounding_box_restriction!(
        x, y, emb_graph::EmbeddedGraph, anchors::AbstractVector{Int}, box_margin::Int
    )
    center = _graph_center_vertex(emb_graph)
    model_coords = _model_frame_coords(emb_graph, center, anchors)
    x_lower, x_upper, y_lower, y_upper = _tightened_bounding_box(
        model_coords, box_margin
    )

    proposed_bounds = Tuple{VariableRef, Float64, Float64}[]
    conflicts = BoundingBoxConflict[]
    for vertex in 1:nv(emb_graph.graph)
        for (var, lower, upper) in (
                (x[vertex], x_lower, x_upper),
                (y[vertex], y_lower, y_upper),
            )
            proposed = _proposed_variable_bounds(var, lower, upper)
            push!(proposed_bounds, (
                var,
                proposed.proposed_lower,
                proposed.proposed_upper,
            ))
            if proposed.proposed_lower > proposed.proposed_upper
                push!(conflicts, (
                    variable = name(var),
                    original_lower = proposed.original_lower,
                    original_upper = proposed.original_upper,
                    proposed_lower = proposed.proposed_lower,
                    proposed_upper = proposed.proposed_upper,
                ))
            end
        end
    end

    isempty(conflicts) || throw(BoundingBoxCandidateRejected(conflicts))
    for (var, lower, upper) in proposed_bounds
        set_lower_bound(var, lower)
        set_upper_bound(var, upper)
    end
    return nothing
end


function _normalize_contraction_vertices(contraction_vertices, n::Int)
    contraction_vertices === nothing && return nothing
    contraction_vertices isa Tuple && length(contraction_vertices) == 2 ||
        throw(ArgumentError("contraction_vertices must be a tuple `(u, v)`"))

    u, v = contraction_vertices
    u isa Integer && v isa Integer ||
        throw(ArgumentError("contraction_vertices must contain integer vertex indices"))
    u_int = Int(u)
    v_int = Int(v)
    1 <= u_int <= n && 1 <= v_int <= n || throw(ArgumentError(
        "contraction_vertices must lie in 1:$n, got ($u_int, $v_int)"
    ))
    u_int != v_int || throw(ArgumentError("contraction_vertices must be distinct"))
    return minmax(u_int, v_int)
end


function _has_common_neighbor(graph::Graphs.AbstractGraph, u::Int, v::Int)::Bool
    return any(neighbor -> has_edge(graph, neighbor, v), neighbors(graph, u))
end


function _validate_contraction_pair(
        pair::NTuple{2, Int}, graph::Graphs.AbstractGraph, coords::AbstractVector{IPoint}
    )::NTuple{2, Int}
    u, v = pair
    coords[u] != coords[v] || throw(ArgumentError(
        "contraction_vertices must have distinct coordinates"
    ))
    !_has_common_neighbor(graph, u, v) || throw(ArgumentError(
        "contraction_vertices must not have a common neighbor"
    ))
    return pair
end


function _contraction_candidates(
        graph::Graphs.AbstractGraph, coords::AbstractVector{IPoint};
        require_nonadjacent::Bool
    )::Vector{NTuple{2, Int}}
    n = nv(graph)
    candidates = NTuple{2, Int}[]
    for u in 1:(n - 1)
        for v in (u + 1):n
            coords[u] == coords[v] && continue
            _has_common_neighbor(graph, u, v) && continue
            require_nonadjacent && has_edge(graph, u, v) && continue
            push!(candidates, (u, v))
        end
    end
    return candidates
end


function _select_contraction_vertices(
        rng::AbstractRNG, graph::Graphs.AbstractGraph, coords::AbstractVector{IPoint}
    )::NTuple{2, Int}
    candidates = _contraction_candidates(graph, coords; require_nonadjacent = true)
    isempty(candidates) &&
        (candidates = _contraction_candidates(graph, coords; require_nonadjacent = false))
    isempty(candidates) && throw(ArgumentError(
        "cannot select contraction vertices without two distinct-coordinate vertices"
    ))
    return rand(rng, candidates)
end


function _resolve_contraction_vertices(
        rng::AbstractRNG, emb_graph::EmbeddedGraph, contraction_vertices
    )::NTuple{2, Int}
    explicit_pair = _normalize_contraction_vertices(contraction_vertices, nv(emb_graph.graph))
    explicit_pair === nothing &&
        return _select_contraction_vertices(rng, emb_graph.graph, emb_graph.coords)
    return _validate_contraction_pair(explicit_pair, emb_graph.graph, emb_graph.coords)
end


function _apply_vertex_contraction!(
        model::Model, x, y, pair::NTuple{2, Int}
    )
    u, v = pair
    @constraint(model, x[u] == x[v])
    @constraint(model, y[u] == y[v])
    return nothing
end


"""
    generate_likely_infeasible_embedding_instance(n; strategy=:bounding_box,
        base=:globally_rigid, seed=0, num_anchors=0, alpha=0.0,
        box_margin=1, contraction_vertices=nothing, kwargs...)

Generate a graph embedding model from a feasible sampled embedding and apply a
structural modification intended to make the model likely infeasible.

The bounding-box strategy moves every face of the tight sampled-coordinate box
inward by `box_margin` and intersects that box with the existing variable
domains. It throws [`BoundingBoxCandidateRejected`](@ref) rather than creating
an invalid variable interval when an intersection is empty.

Supported `strategy` values are `:bounding_box` and `:vertex_contraction`
(`:contraction` is accepted as an alias). Supported `base` values are
`:globally_rigid`, `:laman`, and `:random_2_connected` (`:two_connected` is
accepted as an alias). Extra keyword arguments are forwarded to the selected
base graph sampler.
"""
function generate_likely_infeasible_embedding_instance(
        n::Int; strategy::Symbol = :bounding_box, base::Symbol = :globally_rigid,
        seed::Int = 0, num_anchors::Int = 0, alpha::Real = 0.0,
        box_margin::Integer = 1, contraction_vertices = nothing, kwargs...
    )
    normalized_strategy = validate_infeasible_strategy(strategy)
    validate_infeasible_base(base)
    validate_num_anchors(num_anchors, n)
    alpha_float = validate_inexact_alpha(alpha)
    margin_int = validate_box_margin(box_margin)

    rng = rng_from_seed(seed)
    emb_graph = _sample_infeasible_base(rng, n, base; kwargs...)
    anchors = select_anchor_vertices(rng, n, num_anchors)
    model, x, y = build_embedding_model(emb_graph, anchors; alpha = alpha_float)

    if normalized_strategy == :bounding_box
        _apply_bounding_box_restriction!(x, y, emb_graph, anchors, margin_int)
    else
        pair = _resolve_contraction_vertices(rng, emb_graph, contraction_vertices)
        _apply_vertex_contraction!(model, x, y, pair)
    end

    return model, x, y
end
