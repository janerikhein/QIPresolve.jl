"""
    PresolvingCore

Provide the core quadratic-presolve data structures and transformations.

This module collects the quadratic-expression model, parity presolve pipeline,
and display utilities used throughout QIPresolve.

# See also
- [`QuadExpr`](@ref)
- [`QPModel`](@ref)
- [`normalize!`](@ref)
- [`fix_vars!`](@ref)
"""
module PresolvingCore

include("model/quad_expr.jl")
include("model/constraint.jl")
include("residue/interaction_graph.jl")
include("parity/xor_constraint.jl")
include("parity/propagation.jl")
include("parity/parity_postsolve.jl")
include("model/qp_model.jl")
include("parity/xor_model.jl")
include("parity/parity_presolve.jl")

include("utils.jl")

export
    # Model definition
    QuadExpr,
    IntVar,
    Constraint,
    QPModel,
    InteractionGraph,
    InteractionComponent,
    LinSingleton,
    QuadSingleton,
    NonSingleton,
    decompose,
    # Model transformations
    fix_vars!,
    normalize!


end # module
