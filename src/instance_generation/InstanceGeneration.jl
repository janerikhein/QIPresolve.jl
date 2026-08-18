"""
    InstanceGeneration

Provide graph-generation and embedding-model utilities.

This module exposes model generators for Laman, globally rigid, random
2-connected graph embeddings, and likely infeasible embedding variants.
"""
module InstanceGeneration

export
    # Graph embedding instance generation
    generate_laman_instance,
    generate_globally_rigid_instance,
    generate_2_connected_instance,
    generate_likely_infeasible_embedding_instance,
    # Random QIP generation
    generate_random_qip_model


include("embedded_graph.jl")
include("rigid_gen.jl")
include("random_gen.jl")
include("infeasible_gen.jl")
include("random_qip.jl")

end # module
