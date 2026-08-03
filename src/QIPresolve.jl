"""
    QIPresolve

Load the top-level public API for QIPresolve.

Assemble the presolving core, instance-generation helpers, model I/O utilities, and
package configuration into the main package module.

# See also
- [`PresolvingCore`](@ref)
- [`ModelIO`](@ref)
- [`InstanceGeneration`](@ref)
"""
module QIPresolve

include("config.jl")
include("core/PresolvingCore.jl")
include("instance_generation/InstanceGeneration.jl")
include("model_io/ModelIO.jl")

using .ModelIO: load_moi_model, from_moi, build_moi_model, build_model
using .PresolvingCore: affine_transform!, lin_transform!, fix_vars!, normalize!, presolve!, residue_presolve!, postsolve, PresolveResult

export
    # Presolving Core functionality
    fix_vars!,
    normalize!,
    PresolveResult,
    presolve!,
    residue_presolve!,
    postsolve,
    # Model Input/Output
    load_moi_model,
    from_moi,
    build_moi_model,
    build_model


end
