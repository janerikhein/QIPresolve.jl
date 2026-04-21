"""
    ModelIO

Provide model construction and MathOptInterface import helpers.

This module contains the builder utilities and MOI-file parsing logic used to
construct `PresolvingCore.QPModel` instances from external representations.
"""
module ModelIO

export
    load_moi_model,
    from_moi,
    save_moi,
    build_model

include("qp_model_builder.jl")
include("mof_parser.jl")

end # module
