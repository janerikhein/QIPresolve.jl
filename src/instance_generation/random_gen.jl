using Random
using Graphs


"""
    cycle_edge(u::Int, v::Int, n::Int) -> Bool

Return true if `(u, v)` is an edge of the cycle graph `C_n`.
Assumes `1 <= u < v <= n`.
"""
@inline cycle_edge(u::Int, v::Int, n::Int)::Bool = (v == u + 1) || (u == 1 && v == n)


"""
    _random_2_connected_graph(rng, n; R=10, edge_density=0.5, max_coord_tries=10_000)

Internal RNG-sharing implementation for random 2-connected graph generation.

The graph is constructed from the cycle graph `C_n` and then augmented with an
exact number of extra edges determined by `edge_density`.
"""
function _random_2_connected_graph(
        rng::AbstractRNG, n::Int; R::Int = 10, edge_density::Float64 = 0.5,
        max_coord_tries::Int = 10_000
    )
    n ≥ 3 || throw(ArgumentError("n must be ≥ 3"))
    validate_coordinate_box(n, R)
    0.0 ≤ edge_density ≤ 1.0 || throw(ArgumentError("edge_density must satisfy 0.0 ≤ edge_density ≤ 1.0"))
    max_coord_tries ≥ 1 || throw(ArgumentError("max_coord_tries must be ≥ 1"))

    coords = random_distinct_points(rng, n, R; max_coord_tries = max_coord_tries)
    g = cycle_graph(n)

    m_min = n
    m_max = n * (n - 1) ÷ 2
    m_target = round(Int, m_min + edge_density * (m_max - m_min))
    n_extra = m_target - m_min

    candidates = NTuple{2, Int}[]
    sizehint!(candidates, m_max - m_min)
    for u in 1:(n - 1)
        for v in (u + 1):n
            cycle_edge(u, v, n) && continue
            push!(candidates, (u, v))
        end
    end

    if n_extra > 0
        order = randperm(rng, length(candidates))
        @inbounds for idx in order[1:n_extra]
            u, v = candidates[idx]
            add_edge!(g, u, v)
        end
    end

    return g, coords
end


"""
    random_2_connected_graph(n; R=10, edge_density=0.5, seed=0, max_coord_tries=10_000)

Return `(g::Graph, coords::Vector{IPoint})` where `g` is a simple 2-connected
graph on `n` vertices and `coords` are distinct integer points in `[-R, R]^2`.

The graph is constructed from the cycle graph `C_n` and then augmented with an
exact number of extra edges determined by `edge_density`.
"""
function random_2_connected_graph(
        n::Int; R::Int = 10, edge_density::Float64 = 0.5, seed::Int = 0,
        max_coord_tries::Int = 10_000
    )
    rng = rng_from_seed(seed)
    return _random_2_connected_graph(
        rng, n;
        R = R,
        edge_density = edge_density,
        max_coord_tries = max_coord_tries,
    )
end


"""
    plot_2_connected_graph(g::Graph, coords::Vector{IPoint}; show_labels=false, show_coords=false, kwargs...)

Plot a 2-connected graph with fixed integer coordinates using GraphPlot.jl.
Keyword arguments are forwarded to `GraphPlot.gplot`.
"""
function plot_2_connected_graph(g::Graph, coords::Vector{IPoint}; kwargs...)
    return plot_embedded_graph(g, coords; kwargs...)
end


"""
    generate_2_connected_instance(n; R=10, edge_density=0.5, seed=0,
                                  max_coord_tries=10_000, num_anchors=0,
                                  alpha=0.0)

Generate a 2-connected graph instance and return the corresponding embedding
model as `(model, x, y)`. A seeded random subset of `num_anchors` vertices is
fixed at its generated coordinates. When `alpha > 0`, edge distance equalities
are relaxed to relative intervals.
"""
function generate_2_connected_instance(
        n::Int; R::Int = 10, edge_density::Float64 = 0.5, seed::Int = 0,
        max_coord_tries::Int = 10_000, num_anchors::Int = 0,
        alpha::Real = 0.0
    )
    validate_num_anchors(num_anchors, n)
    alpha_float = validate_inexact_alpha(alpha)
    rng = rng_from_seed(seed)
    g, coords = _random_2_connected_graph(
        rng, n;
        R = R,
        edge_density = edge_density,
        max_coord_tries = max_coord_tries,
    )
    anchors = select_anchor_vertices(rng, n, num_anchors)
    return build_embedding_model(to_embedded(g, coords), anchors; alpha = alpha_float)
end
