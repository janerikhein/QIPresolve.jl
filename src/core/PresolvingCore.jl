module PresolvingCore

include("model/quad_expr.jl")
include("model/constraint.jl")
include("model/qp_model.jl")


#include("gauss_jordan_gf2.jl")
include("parity/xor_constraint.jl")
include("parity/propagation.jl")
include("parity/xor_model.jl")
include("parity/presolve_parity.jl")

include("utils.jl")

export
    # Model definition
    QuadExpr,
    IntVar,
    Constraint,
    QPModel


end # module
