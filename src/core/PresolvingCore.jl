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

using ..PresolveConfig:
    DEFAULT_MAX_PARITY_PHASES,
    DEFAULT_PRESOLVE_PARITY_STRATEGY,
    DEFAULT_PRESOLVE_RESIDUE_STRATEGY,
    DEFAULT_PRESOLVE_RESIDUE_THRESHOLD,
    DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD

include("quad_expr.jl")
include("constraint.jl")
include("interaction_graph.jl")
include("tree_decomposition.jl")
include("xor_constraint.jl")
include("parity_strategy.jl")
include("parity_propagation.jl")
include("parity_postsolve.jl")
include("qp_model.jl")
include("binary_implications.jl")
include("normalization.jl")
include("residue_propagation.jl")
include("xor_model.jl")
include("parity_presolve.jl")
include("presolve.jl")

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
    ConditionedResidueSet,
    ResidueDPResult,
    PresolveResult,
    decompose,
    compute_nonlinear_residue_sets,
    residue_presolve_constraint!,
    residue_presolve!,
    presolve!,
    postsolve,
    # Model transformations
    fix_vars!,
    normalize!


end # module
