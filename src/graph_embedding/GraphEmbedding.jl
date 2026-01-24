module GraphEmbedding

export
    # Laman generation
    random_laman_graph,
    plot_laman_graph,
    generate_laman_instance


include("embedded_graph.jl")
include("laman_gen.jl")

end # module
