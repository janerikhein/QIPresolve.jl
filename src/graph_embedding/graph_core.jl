module GraphCore

export IPoint, EmbeddedGraph, to_embedded

using Graphs
using SimpleWeightedGraphs


"""
    IPoint

Integer 2D point.
"""
struct IPoint
    x::Int
    y::Int
end


"""
    EmbeddedGraph

Graph with embedded coordinates and derived center/farthest metrics.
"""
struct EmbeddedGraph
    graph::SimpleWeightedGraph{Int, Int}
    coords::Vector{IPoint}
    center::Int
    farthest::Int
    dist_from_center::Vector{Float64}
end


"""
    squared_dist2(p1::IPoint, p2::IPoint) -> Int

Squared Euclidean distance between two integer points.
"""
@inline squared_dist2(p1::IPoint, p2::IPoint) = (p1.x - p2.x)^2 + (p1.y - p2.y)^2


"""
    to_embedded(graph::Graph, coords::Vector{IPoint}) -> EmbeddedGraph

Build a weighted graph using squared edge lengths and compute center metadata.
"""
function to_embedded(graph::Graph, coords::Vector{IPoint})::EmbeddedGraph
    n = nv(graph)
    length(coords) == n || throw(ArgumentError("coords must have length nv(graph) = $n, got $(length(coords))"))

    m = ne(graph)
    sources = Vector{Int}(undef, m)
    destinations = Vector{Int}(undef, m)
    sq_edge_lengths = Vector{Int}(undef, m)

    # build edges with squared edge length weights
    @inbounds for (i, e) in enumerate(edges(graph))
        u = src(e)
        v = dst(e)
        sources[i] = u
        destinations[i] = v
        sq_edge_lengths[i] = squared_dist2(coords[u], coords[v])
    end

    # create weighted graph from edges
    wg = SimpleWeightedGraph(sources, destinations, sq_edge_lengths)
    nv(wg) < n && add_vertices!(wg, n - nv(wg))

    pw_dist = pw_shortest_paths(wg)
    center = graph_center(pw_dist)

    # distance vector from center
    dists = Vector{Int}(undef, n)
    @inbounds copyto!(dists, @view pw_dist[center, :])

    # farthest node from center
    farthest = argmax(dists)

    return EmbeddedGraph(wg, copy(coords), center, farthest, dists)
end


"""
    pw_shortest_paths(graph::SimpleWeightedGraph) -> Matrix{Float64}

All-pairs shortest paths using a density-based heuristic to choose the solver.
"""
function pw_shortest_paths(graph::SimpleWeightedGraph)::Matrix{Float64}
    n = nv(graph)
    n == 0 && return Matrix{Float64}(undef, 0, 0)
    n == 1 && return reshape(Float64[0.0], 1, 1)

    m = ne(graph)
    ρ = 2m / (n * (n - 1))  # undirected density in [0,1]

    # Heuristic switch:
    # - Floyd–Warshall for small n or fairly dense graphs
    # - Johnson APSP for sparse graphs
    use_floyd = (n ≤ 400) || (ρ ≥ 0.10)

    distmx = sqrt.(weights(graph))  # edge weights transformed for weighted shortest paths
    if use_floyd && !(distmx isa Matrix)
        distmx = Matrix{Float64}(distmx)
    end

    sp = use_floyd ?
        floyd_warshall_shortest_paths(graph, distmx) :
        johnson_shortest_paths(graph, distmx)

    return Matrix{Float64}(sp.dists)
end


"""
    graph_center(distances::Matrix{Float64}) -> Int

Return the index of a graph center (minimizing maximum distance).
"""
function graph_center(distances::Matrix{Float64})::Int
    n, m = size(distances)
    n == m || throw(ArgumentError("distance matrix must be square, got $(size(distances))"))
    n == 0 && throw(ArgumentError("distance matrix must be non-empty"))

    best_center = 1
    best_radius = Inf

    @inbounds for i in 1:n
        r = maximum(@view distances[i, :])
        if r < best_radius
            best_radius = r
            best_center = i
        end
    end

    return best_center
end


end # module
