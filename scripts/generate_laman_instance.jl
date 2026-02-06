#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuMP: backend
using QIPresolve.GraphEmbedding: generate_laman_instance
using QIPresolve.ModelIO: save_moi



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
    println("Usage: julia scripts/generate_laman_instance.jl [n] [R] [pH2] [seed] [out]")
    return println("Defaults: n=5 R=10 pH2=0.5 seed=0 out=laman_instance")
end

args = copy(ARGS)
if length(args) >= 1 && (args[1] == "-h" || args[1] == "--help")
    usage()
    exit(0)
end

n = length(args) >= 1 ? parse_int(args[1], "n") : 5
R = length(args) >= 2 ? parse_int(args[2], "R") : 100
pH2 = length(args) >= 3 ? parse_float(args[3], "pH2") : 0.5
seed = length(args) >= 4 ? parse_int(args[4], "seed") : 0
out = length(args) >= 5 ? args[5] : "laman_instance.mof.json"

model, _, _ = generate_laman_instance(n; R = R, pH2 = pH2, seed = seed)
path = save_moi(backend(model), out)

println("saved model: $path")
