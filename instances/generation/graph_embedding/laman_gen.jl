module RandomLamanGraphs

export random_laman_graph

using Random
using Graphs
using ..GraphCore: IPoint


@inline function area2(a::IPoint, b::IPoint, c::IPoint)::Int
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
end


@inline is_collinear(a::IPoint, b::IPoint, c::IPoint) = area2(a, b, c) == 0


@inline function rand_point(rng::AbstractRNG, R::Int)::IPoint
    IPoint(rand(rng, -R:R), rand(rng, -R:R))
end


function point_used(coords::Vector{IPoint}, p::IPoint)::Bool
    @inbounds for q in coords
        if q.x == p.x && q.y == p.y
            return true
        end
    end
    return false
end


"Initialize K3 with non-collinear integer coordinates."
function init_triangle(rng::AbstractRNG, R::Int)
    g = Graph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 1, 3)

    # Randomly sample 3 points and check colinearity
    while true
        p1 = rand_point(rng, R)
        p2 = rand_point(rng, R)
        p3 = rand_point(rng, R)
        if (p1 != p2) && (p1 != p3) && (p2 != p3) && !is_collinear(p1, p2, p3)
            return g, IPoint[p1, p2, p3]
        end
    end
end


"""
Henneberg I: add new vertex w connected to two existing vertices a,b.
Returns true on success (graph + coords modified).
"""
function henneberg_I!(graph::Graph, coords::Vector{IPoint}, rng::AbstractRNG, R::Int; max_tries::Int=200)::Bool
    n = nv(graph)

    # sample 2 random nodes
    v1, v2 = rand(rng, 1:n, 2)
    while v2 == v1
        v2 = rand(rng, 1:n)
    end

    pa, pb = coords[v1], coords[v2]
    for _ in 1:max_tries
        p = rand_point(rng, R)
        if !point_used(coords, p) && !is_collinear(pa, pb, p)
            add_vertex!(graph)                 
            push!(coords, p)
            w = nv(graph)
            add_edge!(graph, w, v1)
            add_edge!(graph, w, v2)
            return true
        end
    end
    return false
end

"""
Henneberg II: pick an edge (u,v), remove it, add new vertex w connected to u,v,x.
Returns true on success (graph + coords modified).
"""
function henneberg_II!(g::Graph, coords::Vector{IPoint}, rng::AbstractRNG, R::Int; max_tries::Int=300)::Bool
    n = nv(g)
    ne(g) == 0 && return false

    # pick random existing edge
    edges = collect(edges(g))
    edge = rand(rng, edges)
    u, v = src(edge), dst(edge)

    # choose third attachment x != u,v
    x = rand(rng, 1:n)
    while x == u || x == v
        x = rand(rng, 1:n)
    end

    pu, pv, px = coords[u], coords[v], coords[x]

    for _ in 1:max_tries
        p = rand_point(rng, R)
        if point_used(coords, p)
            continue
        end

        # basic degeneracy rejection with the attachment vertices
        if is_collinear(pu, pv, p) || is_collinear(pu, px, p) || is_collinear(pv, px, p)
            continue
        end

        # apply move
        rem_edge!(g, u, v)
        add_vertex!(g)
        push!(coords, p)
        w = nv(g)
        add_edge!(g, w, u)
        add_edge!(g, w, v)
        add_edge!(g, w, x)
        return true
    end
    return false
end

"""
    random_laman_graph(n; R=10, pH2=0.5, seed=0,
                       max_global_tries=10_000, max_tries_H1=200, max_tries_H2=300)

Return `(g::Graph, coords::Vector{IPoint})` where `g` is (typically) a Laman graph
constructed by Henneberg moves and `coords` are integer points in `[-R,R]^2`.
"""
function random_laman_graph(n::Int; R::Int=10, pH2::Float64=0.5, seed::Int=0,
                            max_global_tries::Int=10_000, max_tries_H1::Int=200, max_tries_H2::Int=300)

    n ≥ 3 || throw(ArgumentError("n must be ≥ 3 (triangle base)"))
    rng = seed == 0 ? Random.default_rng() : MersenneTwister(seed)

    g, coords = init_triangle(rng, R)

    tries = 0

    while nv(g) < n
        tries += 1
        tries > max_global_tries && error("Failed to reach n=$n with R=$R. Try increasing R or tries.")

        useH2 = (rand(rng) < pH2) && (ne(g) > 0)
        ok = useH2 ?
            henneberg_II!(g, coords, rng, R; max_tries=max_tries_H2) :
            henneberg_I!(g, coords, rng, R; max_tries=max_tries_H1)

        ok || continue
    end

    # sanity for Laman counts
    # (will hold combinatorially by construction; this is just a check)
    ne(g) == 2nv(g) - 3 || @warn "Edge count not 2n-3; got $(ne(g)) for n=$(nv(g))"

    return g, coords
end


end # module
