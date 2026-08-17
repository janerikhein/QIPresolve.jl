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
    DEFAULT_PRESOLVE_RESIDUE_STRATEGY,
    DEFAULT_PRESOLVE_RESIDUE_THRESHOLD,
    DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD

include("model/quad_expr.jl")
include("model/constraint.jl")
include("residue/interaction_graph.jl")
include("residue/tree_decomposition.jl")
include("parity/xor_constraint.jl")
include("stats/parity_stats.jl")
include("parity/propagation.jl")
include("parity/parity_postsolve.jl")
include("model/qp_model.jl")
include("stats/model_stats.jl")
include("model/normalization.jl")
include("residue/dp_propagation.jl")
include("parity/xor_model.jl")
include("parity/parity_presolve.jl")
include("presolve.jl")

include("utils.jl")

export
    # Model definition
    QuadExpr,
    IntVar,
    Constraint,
    QPModel,
    ModelStats,
    ParityStats,
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
