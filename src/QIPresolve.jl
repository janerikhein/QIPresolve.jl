"""
    QIPresolve

Load the top-level public API for QIPresolve.

Assemble the presolving core, graph-embedding helpers, model I/O utilities, and
package configuration into the main package module.

# See also
- [`PresolvingCore`](@ref)
- [`ModelIO`](@ref)
- [`GraphEmbedding`](@ref)
"""
module QIPresolve

include("core/PresolvingCore.jl")
include("graph_embedding/GraphEmbedding.jl")
include("model_io/ModelIO.jl")
include("config.jl")

using .ModelIO: load_moi_model, from_moi, build_moi_model, build_model
using .PresolvingCore: affine_transform!, lin_transform!, fix_vars!, normalize!, residue_presolve!

export
    # Presolving Core functionality
    fix_vars!,
    normalize!,
    residue_presolve!,
    # Model Input/Output
    load_moi_model,
    from_moi,
    build_moi_model,
    build_model


end
