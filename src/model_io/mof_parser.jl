import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF


const MOIConSense = Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo, MOI.Interval}
const MOIVarSense = Union{MOI.Integer, MOI.ZeroOne}
const VarSenseMapping = Dict(MOI.Integer => :int, MOI.ZeroOne => :bin)
const ObjSenseMapping = Dict(
    MOI.MIN_SENSE => :min,
    MOI.MAX_SENSE => :max,
    MOI.FEASIBILITY_SENSE => :feas,
)

function load_moi_model(file_path::AbstractString; format::Symbol = :mof)
    if format == :mof
        ff = FF.FORMAT_MOF
    else
        throw(ArgumentError("invalid format identifier: $format"))
    end

    moi = FF.Model(format = ff)
    MOI.read_from_file(moi, file_path)

    return moi
end

function save_moi(model::MOI.ModelLike, filename::AbstractString)
    MOI.write_to_file(model, filename)
    return nothing
end

function from_moi(moi_model)
    qp_model = QPModelBuilder()

    # register objective
    obj_expr_type = MOI.get(moi_model, MOI.ObjectiveFunctionType())
    obj_expr = MOI.get(moi_model, MOI.ObjectiveFunction{obj_expr_type}())
    obj_sense::MOI.OptimizationSense = MOI.get(moi_model, MOI.ObjectiveSense())
    register_objective!(qp_model, obj_expr, obj_sense)

    # register constraints
    for (expr_type, sense_type) in MOI.get(moi_model, MOI.ListOfConstraintTypesPresent())
        con_indices = MOI.get(moi_model, MOI.ListOfConstraintIndices{expr_type, sense_type}())
        for ci in con_indices
            f = MOI.get(moi_model, MOI.ConstraintFunction(), ci)::expr_type
            s = MOI.get(moi_model, MOI.ConstraintSet(), ci)::sense_type
            register_constraint!(qp_model, f, s)
        end
    end

    return qp_model
end


get_bound(moi_sense::MOI.LessThan{T}) where {T <: Real} = (-Inf, Float64(moi_sense.upper))
get_bound(moi_sense::MOI.GreaterThan{T}) where {T <: Real} = (Float64(moi_sense.lower), Inf)
get_bound(moi_sense::MOI.EqualTo{T}) where {T <: Real} = (moi_sense.value, moi_sense.value) .|> Float64
get_bound(moi_sense::MOI.Interval{T}) where {T <: Real} = (moi_sense.lower, moi_sense.upper) .|> Float64

get_term(moi_term::MOI.ScalarAffineTerm{T}) where {T <: Real} = (Float64(moi_term.coefficient), Int(moi_term.variable.value))

# note: somehow mof files scale quad diagonal entries by 2 
function get_term(moi_term::MOI.ScalarQuadraticTerm{T}) where {T <: Real} 
    coeff = Float64(moi_term.coefficient)
    var1 = Int(moi_term.variable_1.value)
    var2 = Int(moi_term.variable_2.value)
    if var1 == var2
        coeff = coeff / 2
    end

    return (coeff, var1, var2)
end


function _parse_quad_expr(expr::MOI.ScalarQuadraticFunction{T}) where {T <: Real}
    quad_terms = Vector{Tuple{Float64, Int, Int}}(undef, length(expr.quadratic_terms))
    @inbounds for (i, term) in enumerate(expr.quadratic_terms)
        quad_terms[i] = get_term(term)
    end

    lin_terms = Vector{Tuple{Float64, Int}}(undef, length(expr.affine_terms))
    @inbounds for (i, term) in enumerate(expr.affine_terms)
        lin_terms[i] = get_term(term)
    end

    return quad_terms, lin_terms
end

function _parse_affine_expr(expr::MOI.ScalarAffineFunction{T}) where {T <: Real}
    lin_terms = Vector{Tuple{Float64, Int}}(undef, length(expr.terms))
    @inbounds for (i, term) in enumerate(expr.terms)
        lin_terms[i] = get_term(term)
    end

    return lin_terms
end

function register_constraint!(qp_model::QPModelBuilder, expr::MOI.ScalarAffineFunction{T}, sense::MOIConSense) where {T <: Real}
    lin_terms = _parse_affine_expr(expr)
    lhs, rhs = get_bound(sense)
    constant = expr.constant

    return register_con!(qp_model; lin_expr_terms = lin_terms, constant = constant, lhs = lhs, rhs = rhs)
end


function register_constraint!(qp_model::QPModelBuilder, expr::MOI.ScalarQuadraticFunction{T}, sense::MOIConSense) where {T <: Real}
    quad_terms, lin_terms = _parse_quad_expr(expr)
    lhs, rhs = get_bound(sense)
    constant = expr.constant

    return register_con!(qp_model; lin_expr_terms = lin_terms, quad_expr_terms = quad_terms, constant = constant, lhs = lhs, rhs = rhs)
end


function register_constraint!(qp_model::QPModelBuilder, expr::MOI.VariableIndex, sense::MOIConSense)
    lhs, rhs = get_bound(sense)

    return register_var_info!(qp_model, expr.value; lb = lhs, ub = rhs)
end


register_constraint!(qp_model::QPModelBuilder, expr::MOI.VariableIndex, ::S) where {S <: MOIVarSense} =
    register_var_info!(qp_model, expr.value; var_type = VarSenseMapping[S])


function register_objective!(qp_model::QPModelBuilder, expr::MOI.ScalarAffineFunction{T}, sense::MOI.OptimizationSense) where {T <: Real}
    lin_terms = _parse_affine_expr(expr)
    constant = expr.constant
    obj_sense = ObjSenseMapping[sense]

    return register_obj!(qp_model, constant, obj_sense; lin_expr_terms = lin_terms)
end
