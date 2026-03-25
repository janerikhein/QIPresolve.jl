#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Graphs
include(joinpath(@__DIR__, "..", "src", "core", "parity", "propagation.jl"))

function show_manager(manager::PropagationManager, title::String)
    println(title)
    println("nvars=$(manager.nvars) nreps=$(manager.nreps)")
    println("labels=$(manager.lit_labels)")
    println("scc_pos_to_rep_pos=$(manager.scc_pos_to_rep_pos)")
    println("edges:")
    for edge in edges(manager.graph)
        println("  $(src(edge)) -> $(dst(edge))")
    end
    println()
    return nothing
end

manager = PropagationManager([1, 2, 3, 4])

add_implication!(manager, VarLit(1, true), VarLit(2, true))
add_implication!(manager, VarLit(2, true), VarLit(3, true))
add_implication!(manager, VarLit(3, true), VarLit(1, true))
add_implication!(manager, VarLit(4, true), VarLit(1, true))
add_implication!(manager, VarLit(3, true), VarLit(4, false))
#add_implication!(manager, VarLit(2, false), VarLit(3, true))
#add_equivalence!(manager, VarLit(3, false), VarLit(4, false))
#fix_var!(manager, 1, true)

show_manager(manager, "Before update!")

update!(manager)

show_manager(manager, "After update!")
