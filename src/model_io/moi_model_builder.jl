using ..PresolvingCore: QPModel, QuadExpr, IntVar, is_binary, vars, get_lin_coeff, get_quad_coeff
import MathOptInterface as MOI
import SCIP

const QPObjSenseMapping = Dict(
    :min => MOI.MIN_SENSE,
    :max => MOI.MAX_SENSE,
    :feas => MOI.FEASIBILITY_SENSE,
    :undef => MOI.FEASIBILITY_SENSE,
)

"""
    build_moi_model(model, solver=nothing)

Build a MathOptInterface model from a `QPModel`.

# Arguments
- `model`: Source `QPModel`.
- `solver`: Optional backend selector. Use `nothing` for an in-memory
  `MOI.Utilities.Model{Float64}` or `:scip` for a SCIP-backed model.

# Returns
- An `MOI.ModelLike` populated with variables, domains, bounds, constraints,
  and objective data from `model`.
"""
function build_moi_model(model::QPModel, solver::Union{Nothing, Symbol} = nothing)
    moi_model = _moi_model_backend(model, solver)
    variable_map = _register_moi_variables!(moi_model, model.vars)

    for con in model.cons
        _register_moi_constraint!(moi_model, variable_map, con.qe, con.lhs, con.rhs)
    end

    if model.infeasible
        infeasible = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0)
        MOI.add_constraint(moi_model, infeasible, MOI.LessThan(-1.0))
    end

    obj_sense = get(QPObjSenseMapping, model.obj_sense, nothing)
    obj_sense === nothing && throw(ArgumentError("invalid objective sense: $(model.obj_sense)"))
    MOI.set(moi_model, MOI.ObjectiveSense(), obj_sense)

    obj_expr = _moi_scalar_function(model.obj_expr, variable_map)
    MOI.set(moi_model, MOI.ObjectiveFunction{typeof(obj_expr)}(), obj_expr)

    return moi_model
end

function _moi_model_backend(::QPModel, ::Nothing)
    return MOI.Utilities.Model{Float64}()
end

function _moi_model_backend(model::QPModel, solver::Symbol)
    if solver == :scip
        _has_quadratic_terms(model.obj_expr) &&
            throw(ArgumentError("SCIP backend does not support quadratic objectives in build_moi_model"))

        return MOI.instantiate(SCIP.Optimizer; with_bridge_type = Float64)
    else
        throw(ArgumentError("unsupported solver: $solver"))
    end
end

function _register_moi_variables!(moi_model::MOI.ModelLike, vars::Dict{VarId, IntVar})
    variable_map = Dict{VarId, MOI.VariableIndex}()
    sizehint!(variable_map, length(vars))

    for var_id in sort!(collect(keys(vars)))
        variable = MOI.add_variable(moi_model)
        variable_map[var_id] = variable
        MOI.set(moi_model, MOI.VariableName(), variable, "x$var_id")

        var_info = vars[var_id]
        if is_binary(var_info)
            MOI.add_constraint(moi_model, variable, MOI.ZeroOne())
        else
            MOI.add_constraint(moi_model, variable, MOI.Integer())
        end
        _register_moi_variable_bounds!(moi_model, variable, var_info.lb, var_info.ub)
    end

    return variable_map
end

function _register_moi_variable_bounds!(moi_model::MOI.ModelLike, variable::MOI.VariableIndex, lb::Float64, ub::Float64)
    bound_set = _moi_bound_set(lb, ub)
    bound_set === nothing && return nothing

    MOI.add_constraint(moi_model, variable, bound_set)
    return nothing
end

function _register_moi_constraint!(
        moi_model::MOI.ModelLike,
        variable_map::Dict{VarId, MOI.VariableIndex},
        expr::QuadExpr,
        lhs::Float64,
        rhs::Float64,
    )
    bound_set = _moi_bound_set(lhs, rhs)
    bound_set === nothing && return nothing

    moi_expr = _moi_scalar_function(expr, variable_map)
    MOI.add_constraint(moi_model, moi_expr, bound_set)
    return nothing
end

function _moi_bound_set(lhs::Float64, rhs::Float64)
    if lhs == -Inf && rhs == Inf
        return nothing
    elseif lhs == rhs
        return MOI.EqualTo(lhs)
    elseif isfinite(lhs) && isfinite(rhs)
        return MOI.Interval(lhs, rhs)
    elseif isfinite(lhs)
        return MOI.GreaterThan(lhs)
    elseif isfinite(rhs)
        return MOI.LessThan(rhs)
    else
        return nothing
    end
end

function _moi_scalar_function(expr::QuadExpr, variable_map::Dict{VarId, MOI.VariableIndex})
    affine_terms = _moi_affine_terms(expr, variable_map)
    quadratic_terms = _moi_quadratic_terms(expr, variable_map)

    if isempty(quadratic_terms)
        return MOI.ScalarAffineFunction(affine_terms, expr.constant)
    else
        return MOI.ScalarQuadraticFunction(quadratic_terms, affine_terms, expr.constant)
    end
end

function _moi_affine_terms(expr::QuadExpr, variable_map::Dict{VarId, MOI.VariableIndex})
    affine_terms = MOI.ScalarAffineTerm{Float64}[]

    for var_id in sort!(collect(vars(expr)))
        coeff = get_lin_coeff(expr, var_id)
        isapprox(coeff, 0.0) && continue
        push!(affine_terms, MOI.ScalarAffineTerm(coeff, _moi_variable(variable_map, var_id)))
    end

    return affine_terms
end

function _moi_quadratic_terms(expr::QuadExpr, variable_map::Dict{VarId, MOI.VariableIndex})
    var_ids = sort!(collect(vars(expr)))
    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]

    for (i, var_id_1) in enumerate(var_ids)
        for var_id_2 in @view var_ids[i:end]
            coeff = get_quad_coeff(expr, var_id_1, var_id_2)
            isapprox(coeff, 0.0) && continue
            moi_coeff = var_id_1 == var_id_2 ? 2.0 * coeff : coeff
            push!(
                quadratic_terms,
                MOI.ScalarQuadraticTerm(
                    moi_coeff,
                    _moi_variable(variable_map, var_id_1),
                    _moi_variable(variable_map, var_id_2),
                ),
            )
        end
    end

    return quadratic_terms
end

function _has_quadratic_terms(expr::QuadExpr)
    var_ids = sort!(collect(vars(expr)))
    for (i, var_id_1) in enumerate(var_ids)
        for var_id_2 in @view var_ids[i:end]
            isapprox(get_quad_coeff(expr, var_id_1, var_id_2), 0.0) || return true
        end
    end

    return false
end

function _moi_variable(variable_map::Dict{VarId, MOI.VariableIndex}, var_id::VarId)
    haskey(variable_map, var_id) || throw(ArgumentError("expression references unknown variable id: $var_id"))
    return variable_map[var_id]
end
