#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Graphs
using Compose
using Cairo
using Fontconfig
using QIPresolve.GraphEmbedding.LamanGeneration: random_laman_graph, plot_laman_graph

function parse_int(arg::String, name::String)::Int
    try
        return parse(Int, arg)
    catch
        error("Invalid $name: $arg")
    end
end

function parse_float(arg::String, name::String)::Float64
    try
        return parse(Float64, arg)
    catch
        error("Invalid $name: $arg")
    end
end

function usage()
    println("Usage: julia scripts/generate_laman_graph.jl [n] [R] [pH2] [seed] [out]")
    println("Defaults: n=5 R=10 pH2=0.5 seed=0 out=laman_graph.png")
end

args = copy(ARGS)
if length(args) >= 1 && (args[1] == "-h" || args[1] == "--help")
    usage()
    exit(0)
end

n = length(args) >= 1 ? parse_int(args[1], "n") : 5
R = length(args) >= 2 ? parse_int(args[2], "R") : 10
pH2 = length(args) >= 3 ? parse_float(args[3], "pH2") : 0.5
seed = length(args) >= 4 ? parse_int(args[4], "seed") : 0
out = length(args) >= 5 ? args[5] : "laman_graph.png"

g, coords = random_laman_graph(n; R=R, pH2=pH2, seed=seed)

println("n=$(nv(g)) m=$(ne(g)) R=$R pH2=$pH2 seed=$seed")
println("coords:")
for (i, p) in enumerate(coords)
    println("$i $(p.x) $(p.y)")
end

println("edges:")
for e in edges(g)
    println("$(src(e)) $(dst(e))")
end

plot_ctx = plot_laman_graph(g, coords; show_coords=true)
draw(PNG(out, 800px, 800px), plot_ctx)
println("saved plot: $out")
