import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF

using .PresolvingCore: QPModel, register_con!, register_var_info!, register_obj!


const MOIConSense = Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo, MOI.Interval}
const MOIVarSense = Union{MOI.Integer, MOI.ZeroOne}
const VarSenseMapping = Dict{MOI.Integer => :int, MOI.ZeroOne => :bin}

function load_model_from_mof(instance_name)
    mof = FF.Model(format = FF.FORMAT_MOF)
    MOI.read_from_file(mof, "./$instance_name.mof.json")
    return mof
end

#TODO: add type casting

_get_bound(moi_sense::MOI.LessThan) = (nothing, moi_sense.upper)
_get_bound(moi_sense::MOI.GreaterThan) = (moi_sense.lower, nothing)
_get_bound(moi_sense::MOI.EqualTo) = (moi_sense.value, moi_sense.value)
_get_bound(moi_sense::MOI.Interval) = (moi_sense.lower, moi_sense.upper)

_get_term(moi_term::MOI.ScalarAffineTerm) = (moi_term.coefficient, moi_term.variable.value)
_get_term(moi_term::MOI.ScalarQuadraticTerm) = (moi_term.coefficient, moi_term.variable_1.value, moi_term.variable_2.value)


function register_constraints!(qp_model::QPModel, moi_model)
    for (expr_type, sense_type) in MOI.get(moi_model, MOI.ListOfConstraintTypesPresent())
        con_indices = MOI.get(moi_model, MOI.ListOfConstraintIndices{expr_type, sense_type}())
        for ci in con_indices
            f = MOI.get(moi_model, MOI.ConstraintFunction(), ci)::expr_type
            s = MOI.get(moi_model, MOI.ConstraintSet(), ci)::sense_type
            _register_constraint!(qp_model, f, s)
        end
    end
end


function _register_constraint!(qp_model::QPModel{T}, expr::E, sense::S) where {E<:MOI.ScalarAffineFunction, S<:MOIConSense}

    lin_terms = Vector{Float64, Int}(undef, length(expr.terms))
    @inbounds for (i, term) in enumerate(expr.terms)
        lin_terms[i] = _get_term(term)
    end

    lhs, rhs = _get_bound(sense)
    constant = expr.constant
    register_con!(qp_model, lin_expr_terms=lin_terms, constant=constant, lhs=lhs, rhs=rhs)
end


function _register_constraint!(qp_model::QPModel, expr::E, sense::S) where {E<:MOI.ScalarQuadraticFunction, S}
    
    quad_terms = Vector{Float64, Int, Int}(undef, length(expr.quadratic_terms))
    @inbounds for (i, term) in enumerate(expr.quadratic_terms)
        quad_terms[i] = _get_term(term)
    end

    lin_terms = Vector{Float64, Int}(undef, length(expr.affine_terms))
    @inbounds for (i, term) in enumerate(expr.affine_terms)
        lin_terms[i] = _get_term(term)
    end

    lhs, rhs = _get_bound(sense)
    constant = expr.constant
    register_con!(qp_model, lin_expr_terms=lin_terms, quad_expr_terms=quad_terms, constant=constant, lhs=lhs, rhs=rhs)
end


function _register_constraint!(qp_model::QPModel, expr::E, sense::S) where {E<:MOI.VariableIndex, S<:MOIConSense}
    lhs, rhs = _get_bound(sense)

    register_var_info!(qp_model, expr.value, lb = lhs, ub = rhs)
end


function _register_constraint!(qp_model::QPModel, expr::E, ::S) where {E<:MOI.VariableIndex, S<:MOIVarSense}
    register_var_info!(qp_model, expr.value, var_type = VarSenseMapping[S])
end




