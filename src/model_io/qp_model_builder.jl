using ..PresolvingCore: QPModel, QuadExpr, Constraint, IntVar

const VarId = Int
const _LinDict = Dict{VarId, Float64}
const _QuadDict = Dict{Tuple{VarId, VarId}, Float64}


struct VarInfo
    lb::Float64
    ub::Float64
    var_type::Symbol

    function VarInfo(lb::Float64, ub::Float64, var_type::Symbol)
        lb > ub && throw(ArgumentError("lb must be <= ub. lb:$lb, ub:$ub"))
        var_type in (:cont, :int, :bin) || throw(ArgumentError("invalid var type: $var_type"))

        return new(lb, ub, var_type)
    end
end


mutable struct QuadExprBuilder
    constant::Float64
    lin::_LinDict
    quad::_QuadDict
end

QuadExprBuilder() = QuadExprBuilder(0, _LinDict(), _QuadDict())
QuadExprBuilder(constant::Float64) = QuadExprBuilder(constant, _LinDict(), _QuadDict())


function add_lin!(expr::QuadExprBuilder, id::VarId, val::Float64)
    return expr.lin[id] = get(expr.lin, id, 0.0) + val
end


function add_quad!(expr::QuadExprBuilder, id1::VarId, id2::VarId, val::Float64)
    id1, id2 = id1 > id2 ? (id2, id1) : (id1, id2)
    expr.quad[(id1, id2)] = get(expr.quad, (id1, id2), 0) + val
end

#TODO: modify this function to work with the current implementation of QuadExpr
function build(builder::QuadExprBuilder)
    quad_terms = Vector{Tuple{Float64, VarId, VarId}}()
    sizehint!(quad_terms, length(builder.quad))
    @inbounds for ((id1, id2), val) in builder.quad
        isapprox(val, 0.0) && continue
        push!(quad_terms, (val, id1, id2))
    end

    lin_terms = Vector{Tuple{Float64, VarId}}()
    sizehint!(lin_terms, length(builder.lin))
    @inbounds for (id, val) in builder.lin
        isapprox(val, 0.0) && continue
        push!(lin_terms, (val, id))
    end

    return QuadExpr(quad_terms, lin_terms; constant = builder.constant)
end


struct ConstraintBuilder
    expr::QuadExprBuilder
    lhs::Float64
    rhs::Float64
end

function build(builder::ConstraintBuilder)
    quad_expr = build(builder.expr)

    return Constraint(quad_expr, builder.lhs, builder.rhs)
end

mutable struct QPModelBuilder
    vars::Dict{VarId, VarInfo}
    cons::Vector{ConstraintBuilder}
    obj_expr::QuadExprBuilder
    obj_sense::Symbol
end


QPModelBuilder() = QPModelBuilder(Dict{VarId, VarInfo}(), Vector{ConstraintBuilder}(), QuadExprBuilder(), :undef)

function register_con!(
    model::QPModelBuilder;
    quad_expr_terms::Vector{Tuple{Float64, VarId, VarId}} = Tuple{Float64, VarId, VarId}[],
    lin_expr_terms::Vector{Tuple{Float64, VarId}} = Tuple{Float64, VarId}[],
    constant::Float64 = 0.0,
    lhs::Float64 = -Inf,
    rhs::Float64 = Inf
)   
    quad_expr = QuadExprBuilder(constant)

    for (coeff, var_id_1, var_id_2) in quad_expr_terms
        register_var_info!(model, var_id_1)
        register_var_info!(model, var_id_2)
        add_quad!(quad_expr, var_id_1, var_id_2, coeff)
    end

    for (coeff, var_id) in lin_expr_terms
        register_var_info!(model, var_id)
        add_lin!(quad_expr, var_id, coeff)
    end

    push!(model.cons, ConstraintBuilder(quad_expr, lhs, rhs))
end

function register_var_info!(model::QPModelBuilder, id::VarId; 
    var_type::Symbol = :cont, lb::Float64=-Inf, ub::Float64=Inf
)
    if !haskey(model.vars, id)
        model.vars[id] = VarInfo(lb, ub, var_type)
        return
    else
        var_info = model.vars[id]
    end

    # restrict var_type
    if (var_info.var_type == :bin && (var_type == :int || var_type == :cont)) || (var_info.var_type == :int && var_type == :cont)
        var_type = var_info.var_type
    end

    # restrict bounds
    lb = max(var_info.lb, lb)
    ub = min(var_info.ub, ub)

    model.vars[id] = VarInfo(lb, ub, var_type)
end

function register_obj!(
    model::QPModelBuilder,
    constant::Float64,
    obj_sense::Symbol;
    quad_expr_terms::Vector{Tuple{Float64, VarId, VarId}} = Tuple{Float64, VarId, VarId}[],
    lin_expr_terms::Vector{Tuple{Float64, VarId}} = Tuple{Float64, VarId}[],
)
    quad_expr = QuadExprBuilder(constant)

    for (coeff, var_id_1, var_id_2) in quad_expr_terms
        register_var_info!(model, var_id_1)
        register_var_info!(model, var_id_2)
        add_quad!(quad_expr, var_id_1, var_id_2, coeff)
    end

    for (coeff, var_id) in lin_expr_terms
        register_var_info!(model, var_id)
        add_lin!(quad_expr, var_id, coeff)
    end

    model.obj_expr = quad_expr
    model.obj_sense = obj_sense
end


function build_model(builder::QPModelBuilder)
    
    vars = Dict{VarId, IntVar}()
    sizehint!(vars, length(builder.vars))
    for (var_id, var_info) in builder.vars
        var_info.var_type == :cont && error("continous variables are unsupported")
        vars[var_id] = IntVar(var_info.lb, var_info.ub)
    end
    
    cons = Vector{Constraint}(undef, length(builder.cons))
    @inbounds for (i, con_builder) in enumerate(builder.cons)
        cons[i] = build(con_builder)
    end

    obj_expr = build(builder.obj_expr)
    obj_sense = builder.obj_sense

    return QPModel(vars, cons, obj_expr, obj_sense)
end
