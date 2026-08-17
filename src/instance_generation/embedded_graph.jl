using Graphs
using GraphPlot
using JuMP
using Random
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

Graph with integer coordinates and squared Euclidean edge weights.
"""
struct EmbeddedGraph
    graph::SimpleWeightedGraph{Int, Int}
    coords::Vector{IPoint}
end


"""Return the default RNG for `seed == 0`, otherwise a seeded Mersenne Twister."""
@inline rng_from_seed(seed::Int) = seed == 0 ? Random.default_rng() : MersenneTwister(seed)


"""
    squared_dist2(p1::IPoint, p2::IPoint) -> Int

Squared Euclidean distance between two integer points.
"""
@inline squared_dist2(p1::IPoint, p2::IPoint) = (p1.x - p2.x)^2 + (p1.y - p2.y)^2


"""
    area2(a::IPoint, b::IPoint, c::IPoint) -> Int

Twice the signed area of triangle `(a, b, c)`.
"""
@inline function area2(a::IPoint, b::IPoint, c::IPoint)::Int
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
end


"""
    is_collinear(a::IPoint, b::IPoint, c::IPoint) -> Bool

Return true if the three points are collinear.
"""
@inline is_collinear(a::IPoint, b::IPoint, c::IPoint) = area2(a, b, c) == 0


"""
    rand_point(rng::AbstractRNG, R::Int) -> IPoint

Sample a random integer point in the box `[-R, R]^2`.
"""
@inline function rand_point(rng::AbstractRNG, R::Int)::IPoint
    return IPoint(rand(rng, -R:R), rand(rng, -R:R))
end


"""
    point_used(coords::Vector{IPoint}, p::IPoint) -> Bool

Return true if `p` already appears in `coords`.
"""
function point_used(coords::Vector{IPoint}, p::IPoint)::Bool
    return any(q -> q == p, coords)
end


"""
    random_distinct_points(rng, n, R; max_coord_tries=10_000)

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
    plot_embedded_graph(g, coords; show_labels=false, show_coords=false, kwargs...)

Plot a graph using its fixed integer coordinates.
"""
function plot_embedded_graph(
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
    validate_coordinate_box(n, R)

Validate that `[-R, R]^2` contains enough integer points for `n` vertices.
"""
function validate_coordinate_box(n::Int, R::Int)
    R >= 0 || throw(ArgumentError("R must be >= 0"))
    side = 2R + 1
    side * side >= n || throw(ArgumentError(
        "box [-R, R]^2 has only $(side * side) integer points, fewer than n=$n"
    ))
    return nothing
end


"""
    validate_num_anchors(num_anchors, n)

Validate the number of fixed-coordinate anchor vertices.
"""
function validate_num_anchors(num_anchors::Int, n::Int)
    0 <= num_anchors <= n || throw(ArgumentError(
        "num_anchors must satisfy 0 <= num_anchors <= n=$n, got $num_anchors"
    ))
    return nothing
end


"""
    validate_inexact_alpha(alpha) -> Float64

Validate the relative distance tolerance used by inexact embedding models.
"""
function validate_inexact_alpha(alpha::Real)::Float64
    alpha_float = Float64(alpha)
    isfinite(alpha_float) || throw(ArgumentError("alpha must be finite, got $alpha"))
    alpha_float >= 0.0 || throw(ArgumentError("alpha must be >= 0, got $alpha"))
    return alpha_float
end


"""
    select_anchor_vertices(rng, n, num_anchors) -> Vector{Int}

Select a reproducible random subset of distinct anchor vertices.
"""
function select_anchor_vertices(
        rng::AbstractRNG, n::Int, num_anchors::Int
    )::Vector{Int}
    validate_num_anchors(num_anchors, n)
    num_anchors == 0 && return Int[]
    return randperm(rng, n)[1:num_anchors]
end


"""
    to_embedded(graph::Graph, coords::Vector{IPoint}) -> EmbeddedGraph

Build a weighted graph whose weights are squared Euclidean edge lengths.
"""
function to_embedded(graph::Graph, coords::Vector{IPoint})::EmbeddedGraph
    n = nv(graph)
    length(coords) == n || throw(ArgumentError("coords must have length nv(graph) = $n, got $(length(coords))"))

    m = ne(graph)
    sources = Vector{Int}(undef, m)
    destinations = Vector{Int}(undef, m)
    sq_edge_lengths = Vector{Int}(undef, m)

    @inbounds for (i, e) in enumerate(edges(graph))
        u = src(e)
        v = dst(e)
        sources[i] = u
        destinations[i] = v
        sq_edge_lengths[i] = squared_dist2(coords[u], coords[v])
    end

    weighted_graph = SimpleWeightedGraph(sources, destinations, sq_edge_lengths)
    nv(weighted_graph) < n && add_vertices!(weighted_graph, n - nv(weighted_graph))

    return EmbeddedGraph(weighted_graph, copy(coords))
end


"""
    get_model_params(emb_graph)

Extract vertex indices, edge indices, and squared edge lengths.
"""
function get_model_params(emb_graph::EmbeddedGraph)
    graph = emb_graph.graph
    vertices = 1:nv(graph)
    edge_indices = Vector{NTuple{2, Int}}(undef, ne(graph))
    graph_weights = weights(graph)
    squared_lengths = Dict{NTuple{2, Int}, eltype(graph_weights)}()
    sizehint!(squared_lengths, ne(graph))

    @inbounds for (idx, edge) in enumerate(edges(graph))
        u = src(edge)
        v = dst(edge)
        edge_index = (u, v)
        edge_indices[idx] = edge_index
        squared_lengths[edge_index] = graph_weights[u, v]
    end

    return vertices, edge_indices, squared_lengths
end


"""
    pw_shortest_paths(graph::SimpleWeightedGraph) -> Matrix{Float64}

Compute all-pairs shortest paths using Euclidean edge lengths.
"""
function pw_shortest_paths(graph::SimpleWeightedGraph)::Matrix{Float64}
    n = nv(graph)
    n == 0 && return Matrix{Float64}(undef, 0, 0)
    n == 1 && return reshape(Float64[0.0], 1, 1)

    m = ne(graph)
    density = 2m / (n * (n - 1))
    use_floyd = (n <= 400) || (density >= 0.1)

    distance_matrix = sqrt.(weights(graph))
    if use_floyd && !(distance_matrix isa Matrix)
        distance_matrix = Matrix{Float64}(distance_matrix)
    end

    shortest_paths = use_floyd ?
        floyd_warshall_shortest_paths(graph, distance_matrix) :
        johnson_shortest_paths(graph, distance_matrix)

    return Matrix{Float64}(shortest_paths.dists)
end


"""
    shortest_paths_from(graph, sources) -> Matrix{Float64}

Compute shortest-path distances from every source using Euclidean edge lengths.
Rows follow the order of `sources`.
"""
function shortest_paths_from(
        graph::SimpleWeightedGraph, sources::AbstractVector{Int}
    )::Matrix{Float64}
    n = nv(graph)
    distances = Matrix{Float64}(undef, length(sources), n)
    distance_matrix = sqrt.(weights(graph))

    for (row, source) in enumerate(sources)
        shortest_paths = dijkstra_shortest_paths(graph, source, distance_matrix)
        @inbounds copyto!(@view(distances[row, :]), shortest_paths.dists)
    end

    return distances
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
        radius = maximum(@view distances[i, :])
        if radius < best_radius
            best_radius = radius
            best_center = i
        end
    end

    return best_center
end


"""
    embedding_references(emb_graph, anchors)

Return the effective reference vertices, their fixed coordinates, and their
shortest-path distances. With no anchors, the graph center at `(0, 0)` is used.
"""
function embedding_references(
        emb_graph::EmbeddedGraph, anchors::AbstractVector{Int}
    )
    n = nv(emb_graph.graph)
    all(anchor -> 1 <= anchor <= n, anchors) || throw(ArgumentError("anchor vertex indices must lie in 1:$n"))
    allunique(anchors) || throw(ArgumentError("anchor vertex indices must be distinct"))

    if isempty(anchors)
        pairwise_distances = pw_shortest_paths(emb_graph.graph)
        center = graph_center(pairwise_distances)
        isfinite(maximum(@view pairwise_distances[center, :])) ||
            throw(ArgumentError("embedding graph must be connected"))
        distances = reshape(copy(@view(pairwise_distances[center, :])), 1, n)
        return Int[center], IPoint[IPoint(0, 0)], distances
    end

    reference_vertices = collect(Int, anchors)
    distances = shortest_paths_from(emb_graph.graph, reference_vertices)
    all(isfinite, distances) || throw(ArgumentError("embedding graph must be connected"))
    reference_coords = emb_graph.coords[reference_vertices]
    return reference_vertices, reference_coords, distances
end


"""
    merged_coordinate_bounds(reference_coords, distances)

Intersect the integer coordinate bounds induced by every reference point.
"""
function merged_coordinate_bounds(
        reference_coords::AbstractVector{IPoint}, distances::Matrix{Float64}
    )
    length(reference_coords) == size(distances, 1) || throw(ArgumentError(
        "reference coordinate count must equal the number of distance rows"
    ))

    n = size(distances, 2)
    x_lower = fill(typemin(Int), n)
    x_upper = fill(typemax(Int), n)
    y_lower = fill(typemin(Int), n)
    y_upper = fill(typemax(Int), n)

    @inbounds for (row, point) in enumerate(reference_coords)
        for vertex in 1:n
            distance = distances[row, vertex]
            isfinite(distance) || throw(ArgumentError("all reference distances must be finite"))
            radius = floor(Int, distance)
            x_lower[vertex] = max(x_lower[vertex], point.x - radius)
            x_upper[vertex] = min(x_upper[vertex], point.x + radius)
            y_lower[vertex] = max(y_lower[vertex], point.y - radius)
            y_upper[vertex] = min(y_upper[vertex], point.y + radius)
        end
    end

    all(x_lower .<= x_upper) || error("anchor-derived x bounds are inconsistent")
    all(y_lower .<= y_upper) || error("anchor-derived y bounds are inconsistent")
    return x_lower, x_upper, y_lower, y_upper
end


"""
    build_embedding_model(emb_graph, anchors=Int[]; alpha=0.0)

Build a JuMP model for an integer graph embedding. Anchor vertices are fixed at
their generated coordinates. With zero anchors, the graph-center symmetry
breaking used by the original generator is retained.

When `alpha > 0`, each edge distance equality `f(x) == c` is relaxed to the
bounded inequality `(1 - alpha)c <= f(x) <= (1 + alpha)c`.
"""
function build_embedding_model(
        emb_graph::EmbeddedGraph, anchors::AbstractVector{Int} = Int[];
        alpha::Real = 0.0
    )
    alpha_float = validate_inexact_alpha(alpha)
    vertices, edge_indices, squared_lengths = get_model_params(emb_graph)
    reference_vertices, reference_coords, distances = embedding_references(emb_graph, anchors)
    bound_distances = alpha_float == 0.0 ? distances : sqrt(1.0 + alpha_float) .* distances
    x_lower, x_upper, y_lower, y_upper = merged_coordinate_bounds(reference_coords, bound_distances)
    model = Model()

    @variable(model, x[i in vertices], Int, lower_bound = x_lower[i], upper_bound = x_upper[i])
    @variable(model, y[i in vertices], Int, lower_bound = y_lower[i], upper_bound = y_upper[i])
    @objective(model, Min, 0)

    for (i, j) in edge_indices
        squared_length = Float64(squared_lengths[(i, j)])
        distance_expr = (x[i] - x[j])^2 + (y[i] - y[j])^2
        if alpha_float == 0.0
            @constraint(model, distance_expr == squared_length)
        else
            @constraint(
                model,
                (1.0 - alpha_float) * squared_length <= distance_expr <=
                    (1.0 + alpha_float) * squared_length
            )
        end
    end

    for (row, reference) in enumerate(reference_vertices)
        point = reference_coords[row]
        for vertex in vertices
            vertex == reference && continue
            radius_squared = alpha_float == 0.0 ?
                round(Int, distances[row, vertex]^2) :
                bound_distances[row, vertex]^2
            @constraint(
                model,
                (x[vertex] - point.x)^2 + (y[vertex] - point.y)^2 <= radius_squared
            )
        end
    end

    for (reference, point) in zip(reference_vertices, reference_coords)
        @constraint(model, x[reference] == point.x)
        @constraint(model, y[reference] == point.y)
    end

    if length(reference_vertices) == 1
        reference = only(reference_vertices)
        point = only(reference_coords)
        farthest = argmax(@view distances[1, :])
        @constraint(model, (y[farthest] - point.y) - (x[farthest] - point.x) <= 0)
        @constraint(model, x[farthest] - point.x >= 0)
        @constraint(model, y[farthest] - point.y >= 0)
    end

    return model, x, y
end
