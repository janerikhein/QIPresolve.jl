module QIPresolve

include("core/PresolvingCore.jl")
include("graph_embedding/GraphEmbedding.jl")
include("model_io/ModelIO.jl")
include("config.jl")

using .ModelIO: load_moi_model, from_moi, build_model
using .PresolvingCore: affine_transform!, lin_transform!, fix_vars!, get_parity_model

export
    # Presolving Core functionality
    fix_vars!,
    get_parity_model,
    # Model Input/Output
    load_moi_model,
    from_moi,
    build_model


end
