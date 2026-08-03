"""
    InstanceGeneration

Provide graph-generation and embedding-model utilities.

This module exposes helpers for generating Laman and 2-connected graph
instances together with JuMP models for their integer-coordinate embeddings.
"""
module InstanceGeneration

export
    # Laman generation
    random_laman_graph,
    plot_laman_graph,
    generate_laman_instance,
    # Random 2-connected generation
    random_2_connected_graph,
    plot_2_connected_graph,
    generate_2_connected_instance,
    # Random QIP generation
    generate_random_qip_model


include("embedded_graph.jl")
include("laman_gen.jl")
include("random_gen.jl")
include("random_qip.jl")

end # module
