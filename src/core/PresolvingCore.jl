module PresolvingCore

include("quad_expr.jl")
include("constraint.jl")
include("model.jl")

export 
    # Model definition
    QuadExpr,
    IntVar,
    Constraint,
    QPModel


end # module
