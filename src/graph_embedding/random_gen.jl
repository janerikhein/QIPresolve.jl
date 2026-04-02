using Random
using Graphs
using GraphPlot


"""
    random_distinct_points(rng::AbstractRNG, n::Int, R::Int; max_coord_tries::Int=10_000) -> Vector{IPoint}

Sample `n` distinct integer points in the box `[-R, R]^2`.
"""
function random_distinct_points(
        rng::AbstractRNG, n::Int, R::Int; max_coord_tries::Int = 10_000
    )::Vector{IPoint}
    coords = Vector{IPoint}()
    sizehint!(coords, n)

    while length(coords) < n
        placed = false
        for _ in 1:max_coord_tries
            p = rand_point(rng, R)
            if !point_used(coords, p)
                push!(coords, p)
                placed = true
                break
            end
        end
        placed || error("Failed to sample $n unique points with R=$R. Try increasing R or max_coord_tries.")
    end

    return coords
end


"""
    cycle_edge(u::Int, v::Int, n::Int) -> Bool

Return true if `(u, v)` is an edge of the cycle graph `C_n`.
Assumes `1 <= u < v <= n`.
"""
@inline cycle_edge(u::Int, v::Int, n::Int)::Bool = (v == u + 1) || (u == 1 && v == n)


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
    n ≥ 3 || throw(ArgumentError("n must be ≥ 3"))
    R ≥ 0 || throw(ArgumentError("R must be ≥ 0"))
    0.0 ≤ edge_density ≤ 1.0 || throw(ArgumentError("edge_density must satisfy 0.0 ≤ edge_density ≤ 1.0"))
    max_coord_tries ≥ 1 || throw(ArgumentError("max_coord_tries must be ≥ 1"))

    side = 2R + 1
    side * side ≥ n || throw(ArgumentError("box [-R, R]^2 has only $(side * side) integer points, fewer than n=$n"))

    rng = seed == 0 ? Random.default_rng() : MersenneTwister(seed)

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
    plot_2_connected_graph(g::Graph, coords::Vector{IPoint}; show_labels=false, show_coords=false, kwargs...)

Plot a 2-connected graph with fixed integer coordinates using GraphPlot.jl.
Keyword arguments are forwarded to `GraphPlot.gplot`.
"""
function plot_2_connected_graph(
        g::Graph, coords::Vector{IPoint};
        show_labels::Bool = false,
        show_coords::Bool = false,
        kwargs...
    )
    n = nv(g)
    length(coords) == n || throw(ArgumentError("coords must have length nv(g) = $n, got $(length(coords))"))

    xs = [p.x for p in coords]
    ys = [p.y for p in coords]
    labels = show_coords ? [string(p.x, ",", p.y) for p in coords] :
        (show_labels ? collect(1:n) : nothing)
    plot_kwargs = merge((; nodesize = 0.03, nodelabel = labels), (; kwargs...))

    return GraphPlot.gplot(g, xs, ys; plot_kwargs...)
end


"""
    generate_2_connected_instance(n; R=10, edge_density=0.5, seed=0, max_coord_tries=10_000)

Generate a 2-connected graph instance and return the corresponding embedding
model as `(model, x, y)`.
"""
function generate_2_connected_instance(
        n::Int; R::Int = 10, edge_density::Float64 = 0.5, seed::Int = 0,
        max_coord_tries::Int = 10_000
    )
    g, coords = random_2_connected_graph(
        n;
        R = R,
        edge_density = edge_density,
        seed = seed,
        max_coord_tries = max_coord_tries,
    )
    emb_graph = to_embedded(g, coords)
    return build_embedding_model(emb_graph)
end
