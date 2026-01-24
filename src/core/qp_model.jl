const VarId = Int
const _LinDict = Dict{VarId, Float64}
const _QuadDict = Dict{VarId, Dict{VarId, Float64}}


#TODO: make QPModel{T} add copy function to cast to different type when apllicable
mutable struct QPModel
    vars::Dict{VarId, ModelVar}
    cons::Vector{Constraint}
    obj::Objective
end

function QPModel(nvars::Int)
    vars = Dict{VarId, ModelVar}()
    sizehint!(vars, n)
    cons = Vector{Constraint}()
    obj 
end

function register_con!(
    model::QPModel;
    quad_expr_terms::Vector{Tuple{T, VarId, VarId}} = _QuadDict(),
    lin_expr_terms::Vector{Tuple{T, VarId}} = _LinDict(),
    constant::T = zero(T),
    lhs::Union{T, nothing} = nothing,
    rhs::Union{T, nothing} = nothing
) where {T <: Real}
    # add constraint in some nice form
end

function register_var_info!(model::QPModel, id::VarId; 
    var_type::Symbol = :cont, lb::Union{T, nothing}=nothing, ub::Union{T, nothing}=nothing
) where {T<:Real}
    # add var
end

function register_obj!(
    model::QPModel,
    lin_expr_terms::Vector{Tuple{T, VarId}},
    quad_expr_terms::Vector{Tuple{T, VarId, VarId}},
    constant::T;
    obj_sense::Symbol
) where {T <: Real}

end

struct ModelVar
    id::VarId
    lb::Float64
    ub::Float64
    is_int::Bool

    function ModelVar(id::VarId, lb::Float64, ub::Float64, is_int::Bool)
        if lb > ub
            throw(ArgumentError("lb must be <= ub. lb:$lb, ub:$ub"))
        end
        if is_int
            lb = ceil(lb)
            ub = floor(ub)
            if lb > ub
                throw(ArgumentError("ceil(lb) must be <= floor(ub) for int var. ceil(lb):$lb, floor(ub):$ub"))
            end
        end
        return new(id, lb, ub, is_integer)
    end
end

ContVar(id::VarId; lb::Float64 = -Inf, ub::Float64 = Inf) = ModelVar(id, lb, ub, false)
IntVar(id::VarId; lb::Float64 = -Inf, ub::Float64 = Inf) = ModelVar(id, lb, ub, true)
BinVar(id::VarId) = ModelVar(id, 0.0, 1.0, true)


function var_type(var::ModelVar)
    var.is_int || return :cont
    var.lb == 0.0 && var.ub == 1.0 && return :bin
    return :int
end


mutable struct QuadExpr
    constant::Float64
    lin::_LinDict
    quad::_QuadDict
end

QuadExpr() = QuadExpr(0, _LinDict(), _QuadDict())

function add_lin!(expr::QuadExpr, id::VarId, val::Float64)
    return expr.lin[id] = get(expr.lin, id, 0.0) + val
end

function set_lin!(expr::QuadExpr, id::VarId, val::Float64)
    return expr.lin[id] = val
end

function add_quad!(expr::QuadExpr, id1::VarId, id2::VarId, val::Float64; _sym = true)
    inner = get!(expr.quad, id1) do
        Dict{VarId, Float64}()
    end
    inner[id2] = get(inner, id2, 0.0) + val
    return _sym && id1 != id2 && add_quad!(expr, id2, id1, _sym = false)
end

function set_quad!(expr::QuadExpr, id1::VarId, id2::VarId, val::Float64; _sym = true)
    inner = get!(expr.quad, id1) do
        Dict{VarId, Float64}()
    end
    inner[id2] = val
    return _sym && id1 != id2 && set_quad!(expr, id2, id1, _sym = false)
end


mutable struct Objective
    expr::QuadExpr
    sense::Symbol
end

Objective() = Objective(QuadExpr(), :undef)


mutable struct Constraint
    expr::QuadExpr
    lhs::Float64
    rhs::Float64
end


function normalize!(con::Constraint)
    # add expr constant to lhs/rhs
    con.lhs += con.expr.constant
    con.rhs += con.expr.constant
    return con.expr.constant = 0.0
end


function normalize!(model::QPModel)
    # identify offset binary vars and shift to (0,1) valued (use affine var transform)

    # linearize binary quad terms

    # normalize all constraints
end

function affine_var_transform!(model::QPModel, id::VarId, factor::Float64, offset::Float64)

end
