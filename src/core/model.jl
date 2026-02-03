const VarId = Int 

"""
    IntVar(lb, ub)

Integer variable bounds for a `QPModel`.
"""
struct IntVar
    lb::Float64
    ub::Float64
end

"""
    QPModel(vars, cons, obj_expr, obj_sense)

Quadratic program model with integer variables, constraints, and an objective.

`vars` maps variable ids to bounds, `cons` stores constraints, and `obj_expr`
is a `QuadExpr` with sense `obj_sense`.
"""
mutable struct QPModel
    vars::Dict{VarId, IntVar}
    cons::Vector{Constraint}
    obj_expr::QuadExpr
    obj_sense::Symbol
    _max_var_id::VarId

    function QPModel(vars::Dict{VarId, IntVar}, cons::Vector{Constraint}, obj_expr::QuadExpr, obj_sense::Symbol)
        max_var_id = max(keys(vars)...)

        new(vars, cons, obj_expr, obj_sense, max_var_id)
    end
end


"""
    add_var!(model, info) -> VarId

Add a new variable to `model` with bounds `info` and return its id.
"""
function add_var!(model::QPModel, info::IntVar)
    new_id = model._max_var_id += 1
    model.vars[new_id] = info 

    return new_id
end

"""
    set_var_bounds!(model, id, lb, ub)

Set bounds for variable `id` to `[lb, ub]`.
"""
@inline set_var_bounds!(model::QPModel, id::VarId, lb::Float64, ub::Float64) = (model.vars[id] = IntVar(lb, ub))


"""
    var_bound_shift!(model, var_id, shift)

Shift variable `var_id` by a constant in all constraints and the objective.
This mutates `model` and normalizes constraints after the update.
"""
function var_bound_shift!(model::QPModel, var_id::VarId, shift::Float64)
    # apply shift to constraints
    for con in model.cons
        affine_transform!(con.qe, var_id, 1.0, shift)
        normalize!(con)
    end
    # apply shift to objective
    affine_transform!(model.obj_expr, var_id, 1.0, shift)

    # register domain shift
    info = model.vars[var_id]
    set_var_bounds!(model, var_id, info.lb - shift, info.ub - shift)
end


"""
    affine_transform!(model, var_id, scale, offset)

Apply the affine substitution `x_var_id := scale * x_var_id + offset`
to all constraints and the objective in `model`.
"""
function affine_transform!(model::QPModel, var_id::VarId, scale::Float64, offset::Float64)
    # apply transformation to objective expression
    for con in model.cons
        affine_transform!(con, var_id, scale, offset)
    end

    # apply transformation to objective expression
    affine_transform!(model.obj_expr, var_id, scale, offset)
end


"""
    lin_transform!(model, var_id, other_id, a, b)

Apply the linear substitution `x_var_id := a * x_var_id + b * x_other_id`
to all constraints and the objective in `model`.
"""
function lin_transform!(model::QPModel, var_id::VarId, other_id::VarId, a::Float64, b::Float64)
    # apply transformation to objective expression
    for con in model.cons
        for con in model.cons
            lin_transform!(con, var_id, other_id, a, b)
        end
    end

    # apply transformation to objective expression
    lin_transform!(model.obj_expr, var_id, other_id, a, b)
end


