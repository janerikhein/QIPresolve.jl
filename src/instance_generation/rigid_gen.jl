using Graphs
using Random


"""
    init_triangle(rng, R; max_tries=10_000) -> Tuple{Graph, Vector{IPoint}}

Initialize `K_3` with distinct, non-collinear integer coordinates.
"""
function init_triangle(
        rng::AbstractRNG, R::Int; max_tries::Int = 10_000
    )
    for _ in 1:max_tries
        coords = random_distinct_points(rng, 3, R; max_coord_tries = max_tries)
        is_collinear(coords...) && continue
        return complete_graph(3), coords
    end
    error("Failed to initialize a non-collinear triangle with R=$R. Try increasing R or max_tries.")
end


"""
    no_three_collinear(coords) -> Bool

Return true when every triple of points in `coords` is non-collinear.
"""
function no_three_collinear(coords::AbstractVector{IPoint})::Bool
    n = length(coords)
    for i in 1:(n - 2)
        for j in (i + 1):(n - 1)
            for k in (j + 1):n
                is_collinear(coords[i], coords[j], coords[k]) && return false
            end
        end
    end
    return true
end


"""
    init_complete_four(rng, R; max_tries=10_000)

Initialize `K_4` with distinct integer coordinates and no collinear triple.
"""
function init_complete_four(
        rng::AbstractRNG, R::Int; max_tries::Int = 10_000
    )
    for _ in 1:max_tries
        coords = random_distinct_points(rng, 4, R; max_coord_tries = max_tries)
        no_three_collinear(coords) || continue
        return complete_graph(4), coords
    end
    error("Failed to initialize a nondegenerate K_4 with R=$R. Try increasing R or max_tries.")
end


"""
    henneberg_I!(graph, coords, rng, R; max_tries=200) -> Bool

Add a new vertex connected to two existing vertices.
"""
function henneberg_I!(
        graph::Graph, coords::Vector{IPoint}, rng::AbstractRNG, R::Int;
        max_tries::Int = 200
    )::Bool
    n = nv(graph)
    v1, v2 = rand(rng, 1:n, 2)
    while v2 == v1
        v2 = rand(rng, 1:n)
    end

    p1, p2 = coords[v1], coords[v2]
    for _ in 1:max_tries
        point = rand_point(rng, R)
        if !point_used(coords, point) && !is_collinear(p1, p2, point)
            add_vertex!(graph)
            push!(coords, point)
            new_vertex = nv(graph)
            add_edge!(graph, new_vertex, v1)
            add_edge!(graph, new_vertex, v2)
            return true
        end
    end
    return false
end


"""
    henneberg_II!(graph, coords, rng, R; max_tries=300) -> Bool

Perform a two-dimensional 1-extension: remove an edge and add a new vertex
adjacent to both endpoints and a distinct third vertex.
"""
function henneberg_II!(
        graph::Graph, coords::Vector{IPoint}, rng::AbstractRNG, R::Int;
        max_tries::Int = 300
    )::Bool
    n = nv(graph)
    ne(graph) == 0 && return false

    edge = rand(rng, collect(edges(graph)))
    u, v = src(edge), dst(edge)
    third = rand(rng, 1:n)
    while third == u || third == v
        third = rand(rng, 1:n)
    end

    pu, pv, pt = coords[u], coords[v], coords[third]
    for _ in 1:max_tries
        point = rand_point(rng, R)
        point_used(coords, point) && continue
        if is_collinear(pu, pv, point) ||
                is_collinear(pu, pt, point) ||
                is_collinear(pv, pt, point)
            continue
        end

        rem_edge!(graph, u, v)
        add_vertex!(graph)
        push!(coords, point)
        new_vertex = nv(graph)
        add_edge!(graph, new_vertex, u)
        add_edge!(graph, new_vertex, v)
        add_edge!(graph, new_vertex, third)
        return true
    end
    return false
end


function validate_laman_parameters(
        n::Int, R::Int, pH2::Float64, max_global_tries::Int,
        max_tries_H1::Int, max_tries_H2::Int
    )
    n >= 3 || throw(ArgumentError("n must be >= 3 (triangle base)"))
    validate_coordinate_box(n, R)
    0.0 <= pH2 <= 1.0 || throw(ArgumentError("pH2 must satisfy 0.0 <= pH2 <= 1.0"))
    max_global_tries >= 1 || throw(ArgumentError("max_global_tries must be >= 1"))
    max_tries_H1 >= 1 || throw(ArgumentError("max_tries_H1 must be >= 1"))
    max_tries_H2 >= 1 || throw(ArgumentError("max_tries_H2 must be >= 1"))
    return nothing
end


function _random_laman_graph(
        rng::AbstractRNG, n::Int; R::Int = 10, pH2::Float64 = 0.5,
        max_global_tries::Int = 10_000, max_tries_H1::Int = 200,
        max_tries_H2::Int = 300
    )
    validate_laman_parameters(n, R, pH2, max_global_tries, max_tries_H1, max_tries_H2)
    graph, coords = init_triangle(rng, R; max_tries = max_global_tries)

    tries = 0
    while nv(graph) < n
        tries += 1
        tries > max_global_tries && error(
            "Failed to reach n=$n with R=$R. Try increasing R or max_global_tries."
        )

        use_h2 = rand(rng) < pH2
        succeeded = use_h2 ?
            henneberg_II!(graph, coords, rng, R; max_tries = max_tries_H2) :
            henneberg_I!(graph, coords, rng, R; max_tries = max_tries_H1)
        succeeded || continue
    end

    ne(graph) == 2nv(graph) - 3 || error(
        "Laman construction invariant failed: got $(ne(graph)) edges for n=$(nv(graph))"
    )
    return graph, coords
end


"""
    random_laman_graph(n; R=10, pH2=0.5, seed=0, ...)

Generate a Laman graph by Henneberg-I and Henneberg-II operations together
with distinct integer coordinates in `[-R, R]^2`.
"""
function random_laman_graph(
        n::Int; R::Int = 10, pH2::Float64 = 0.5, seed::Int = 0,
        max_global_tries::Int = 10_000, max_tries_H1::Int = 200,
        max_tries_H2::Int = 300
    )
    rng = rng_from_seed(seed)
    return _random_laman_graph(
        rng, n;
        R = R,
        pH2 = pH2,
        max_global_tries = max_global_tries,
        max_tries_H1 = max_tries_H1,
        max_tries_H2 = max_tries_H2,
    )
end


"""Plot a Laman graph with its fixed integer coordinates."""
function plot_laman_graph(g::Graph, coords::Vector{IPoint}; kwargs...)
    return plot_embedded_graph(g, coords; kwargs...)
end


"""
    generate_laman_instance(n; num_anchors=0, ...)

Generate a Laman graph embedding model and fix a seeded random subset of
`num_anchors` vertices at their generated coordinates.
"""
function generate_laman_instance(
        n::Int; R::Int = 10, pH2::Float64 = 0.5, seed::Int = 0,
        max_global_tries::Int = 10_000, max_tries_H1::Int = 200,
        max_tries_H2::Int = 300, num_anchors::Int = 0
    )
    validate_num_anchors(num_anchors, n)
    rng = rng_from_seed(seed)
    graph, coords = _random_laman_graph(
        rng, n;
        R = R,
        pH2 = pH2,
        max_global_tries = max_global_tries,
        max_tries_H1 = max_tries_H1,
        max_tries_H2 = max_tries_H2,
    )
    anchors = select_anchor_vertices(rng, n, num_anchors)
    return build_embedding_model(to_embedded(graph, coords), anchors)
end


function validate_globally_rigid_parameters(
        n::Int, R::Int, max_global_tries::Int, max_tries_H2::Int
    )
    n >= 4 || throw(ArgumentError("n must be >= 4 (K_4 base)"))
    validate_coordinate_box(n, R)
    max_global_tries >= 1 || throw(ArgumentError("max_global_tries must be >= 1"))
    max_tries_H2 >= 1 || throw(ArgumentError("max_tries_H2 must be >= 1"))
    return nothing
end


function _random_globally_rigid_graph(
        rng::AbstractRNG, n::Int; R::Int = 10,
        max_global_tries::Int = 10_000, max_tries_H2::Int = 300
    )
    validate_globally_rigid_parameters(n, R, max_global_tries, max_tries_H2)
    graph, coords = init_complete_four(rng, R; max_tries = max_global_tries)

    tries = 0
    while nv(graph) < n
        tries += 1
        tries > max_global_tries && error(
            "Failed to reach n=$n with R=$R. Try increasing R or max_global_tries."
        )
        henneberg_II!(graph, coords, rng, R; max_tries = max_tries_H2) || continue
    end

    ne(graph) == 2nv(graph) - 2 || error(
        "global-rigidity construction invariant failed: got $(ne(graph)) edges for n=$(nv(graph))"
    )
    return graph, coords
end


"""
    random_globally_rigid_graph(n; R=10, seed=0, ...)

Generate a generically globally rigid graph from `K_4` by repeated
Henneberg-II/1-extension operations, with distinct integer coordinates.
"""
function random_globally_rigid_graph(
        n::Int; R::Int = 10, seed::Int = 0,
        max_global_tries::Int = 10_000, max_tries_H2::Int = 300
    )
    rng = rng_from_seed(seed)
    return _random_globally_rigid_graph(
        rng, n;
        R = R,
        max_global_tries = max_global_tries,
        max_tries_H2 = max_tries_H2,
    )
end


"""Plot a globally rigid graph with its fixed integer coordinates."""
function plot_globally_rigid_graph(g::Graph, coords::Vector{IPoint}; kwargs...)
    return plot_embedded_graph(g, coords; kwargs...)
end


"""
    generate_globally_rigid_instance(n; num_anchors=0, ...)

Generate a globally rigid graph embedding model and fix a seeded random subset
of `num_anchors` vertices at their generated coordinates.
"""
function generate_globally_rigid_instance(
        n::Int; R::Int = 10, seed::Int = 0,
        max_global_tries::Int = 10_000, max_tries_H2::Int = 300,
        num_anchors::Int = 0
    )
    validate_num_anchors(num_anchors, n)
    rng = rng_from_seed(seed)
    graph, coords = _random_globally_rigid_graph(
        rng, n;
        R = R,
        max_global_tries = max_global_tries,
        max_tries_H2 = max_tries_H2,
    )
    anchors = select_anchor_vertices(rng, n, num_anchors)
    return build_embedding_model(to_embedded(graph, coords), anchors)
end
