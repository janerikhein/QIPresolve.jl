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
    is_binary(var)

Check whether `var` has binary domain `[0, 1]`.
"""
is_binary(var::IntVar) = (var.lb == 0.0 && var.ub == 1.0)

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
    infeasible::Bool

    function QPModel(vars::Dict{VarId, IntVar}, cons::Vector{Constraint}, obj_expr::QuadExpr, obj_sense::Symbol)
        max_var_id = maximum(keys(vars); init = 0)

        return new(vars, cons, obj_expr, obj_sense, max_var_id, false)
    end
end

"""
    nvars(model)

Count active variables in `model`.
"""
nvars(model::QPModel) = length(model.vars)

const _EMPTY_QUAD_TERMS = Tuple{Float64, VarId, VarId}[]


"""
    add_var!(model, info) -> VarId

Add a new variable to `model` with bounds `info` and return its id.
"""
function add_var!(model::QPModel, info::IntVar)
    new_id = model._max_var_id += 1
    model.vars[new_id] = info

    return new_id
end

@inline next_con_id(model::QPModel) = maximum((con.id for con in model.cons); init = 0) + 1

@inline lin_expr(lin_terms::Vector{Tuple{Float64, VarId}}) = QuadExpr(_EMPTY_QUAD_TERMS, lin_terms)

"""
    add_constraint!(model, qe, lhs, rhs)

Append a new constraint to `model`.

# Returns
- The newly created `Constraint`.

# Side Effects
- Mutates `model.cons`.
"""
function add_constraint!(model::QPModel, qe::QuadExpr, lhs::Float64, rhs::Float64)
    con = Constraint(next_con_id(model), qe, lhs, rhs)
    push!(model.cons, con)
    return con
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
        affine_transform!(con, var_id, 1.0, shift)
    end
    # apply shift to objective
    affine_transform!(model.obj_expr, var_id, 1.0, shift)

    # register domain shift
    info = model.vars[var_id]
    return set_var_bounds!(model, var_id, info.lb - shift, info.ub - shift)
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
    return affine_transform!(model.obj_expr, var_id, scale, offset)
end

"""
    fix_parities!(model, propagator, postsolver=nothing)

Rewrite variables with fixed parity as `x_var_id := 2 * x_var_id + offset`,
where `offset` is `0` for even and `1` for odd.

# Returns
- The number of transformed variables.

# Side Effects
- Mutates `model`.
- May mutate `postsolver`.
"""
fix_parities!(model::QPModel, propagator::PropagationManager) = fix_parities!(model, propagator, nothing)

function fix_parities!(
    model::QPModel,
    propagator::PropagationManager,
    postsolver::Union{Nothing, ParityPostsolver},
)
    nfixed = 0

    for lit in propagator.pos_to_lit
        lit.neg && continue
        haskey(model.vars, lit.vid) || continue

        val = fixed_value(propagator, lit.vid)
        val === nothing && continue

        offset = val ? 1.0 : 0.0
        var = model.vars[lit.vid]
        new_lb = ceil((var.lb - offset) / 2.0)
        new_ub = floor((var.ub - offset) / 2.0)
        new_lb == 0.0 && (new_lb = 0.0)
        new_ub == 0.0 && (new_ub = 0.0)

        postsolver !== nothing && append_fixed_bit!(postsolver, lit.vid, val)
        affine_transform!(model, lit.vid, 2.0, offset)
        set_var_bounds!(model, lit.vid, new_lb, new_ub)
        nfixed += 1
    end

    return nfixed
end


"""
    lin_transform!(model, var_id, other_id, a, b; c = 0.0)

Apply the linear substitution `x_var_id := a * x_var_id + b * x_other_id + c`
to all constraints and the objective in `model`.
"""
function lin_transform!(
    model::QPModel,
    var_id::VarId,
    other_id::VarId,
    a::Float64,
    b::Float64;
    c::Float64 = 0.0,
)
    # apply transformation to objective expression
    for con in model.cons
        lin_transform!(con, var_id, other_id, a, b; c = c)
    end

    if b != 0.0 && var_id != other_id && has_var(model.obj_expr, var_id) && !has_var(model.obj_expr, other_id)
        add_var!(model.obj_expr, other_id; clear_buf = true)
    end

    # apply transformation to objective expression
    return lin_transform!(model.obj_expr, var_id, other_id, a, b; c = c)
end

fix_parity_patterns!(model::QPModel, propagator::ParityPropagator) = fix_parity_patterns!(model, propagator, nothing)

"""
    fix_parity_patterns!(model, propagator, postsolver=nothing)

Rewrite eligible parity-equivalent SCCs using a fresh binary helper variable.

For each unlabeled positive/negative SCC pair with enough structure, introduce a
binary helper variable, rewrite member variables as affine functions of that
helper, add linking constraints, and update the propagator SCC structure.

# Returns
- The number of rewritten variables.

# Side Effects
- Mutates `model` and `propagator`.
- May mutate `postsolver`.
"""
function fix_parity_patterns!(
    model::QPModel,
    propagator::ParityPropagator,
    postsolver::Union{Nothing, ParityPostsolver},
)
    candidate_vids = VarId[]

    for scc_pos in eachindex(propagator.scc_pos_to_rep_pos)
        propagator.lit_labels[scc_pos] == UNDEF || continue
        rep_pos = propagator.scc_pos_to_rep_pos[scc_pos]
        rep_lit = propagator.pos_to_lit[rep_pos]
        rep_lit.neg && continue

        neg_rep = repr!(propagator, propagator.lit_to_pos[negated(rep_lit)])
        neg_scc = propagator.rep_pos_to_scc_pos[neg_rep]
        propagator.lit_labels[neg_scc] == UNDEF || continue
        _component_size(propagator, rep_pos) >= 2 || continue
        _component_size(propagator, neg_rep) >= 2 || continue

        push!(candidate_vids, rep_lit.vid)
    end

    nrewritten = 0

    for scc_vid in candidate_vids
        haskey(propagator.lit_to_pos, VarLit(scc_vid, false)) || continue

        rep_pos = repr!(propagator, propagator.lit_to_pos[VarLit(scc_vid, false)])
        rep_lit = propagator.pos_to_lit[rep_pos]
        rep_lit.neg && continue

        scc_pos = propagator.rep_pos_to_scc_pos[rep_pos]
        neg_rep = repr!(propagator, propagator.lit_to_pos[negated(rep_lit)])
        neg_scc = propagator.rep_pos_to_scc_pos[neg_rep]
        propagator.lit_labels[scc_pos] == UNDEF || continue
        propagator.lit_labels[neg_scc] == UNDEF || continue

        component_positions = [
            pos for pos in eachindex(propagator.pos_to_lit) if repr!(propagator, pos) == rep_pos
        ]
        length(component_positions) >= 2 || continue

        component_lits = [propagator.pos_to_lit[pos] for pos in component_positions]
        all(haskey(model.vars, lit.vid) for lit in component_lits) || continue
        all(isfinite(model.vars[lit.vid].lb) && isfinite(model.vars[lit.vid].ub) for lit in component_lits) || continue
        @assert length(unique(lit.vid for lit in component_lits)) == length(component_lits)

        b_scc = add_var!(model, IntVar(0.0, 1.0))
        postsolver !== nothing && ensure_tracked_var!(postsolver, b_scc)

        for lit in component_lits
            var = model.vars[lit.vid]
            e_lo = ceil(var.lb / 2.0)
            e_hi = floor(var.ub / 2.0)
            o_lo = ceil((var.lb - 1.0) / 2.0)
            o_hi = floor((var.ub - 1.0) / 2.0)

            if lit.neg
                l0, u0 = o_lo, o_hi
                l1, u1 = e_lo, e_hi
                postsolver !== nothing && append_binary_bit!(postsolver, lit.vid, b_scc; negated = true)
                lin_transform!(model, lit.vid, b_scc, 2.0, -1.0; c = 1.0)
            else
                l0, u0 = e_lo, e_hi
                l1, u1 = o_lo, o_hi
                postsolver !== nothing && append_binary_bit!(postsolver, lit.vid, b_scc)
                lin_transform!(model, lit.vid, b_scc, 2.0, 1.0)
            end

            m = max(0.0, u0 - l1, u1 - l0)
            new_lb = min(l0, l1)
            new_ub = max(u0, u1)
            new_lb == 0.0 && (new_lb = 0.0)
            new_ub == 0.0 && (new_ub = 0.0)

            add_constraint!(
                model,
                lin_expr(Tuple{Float64, VarId}[(1.0, lit.vid), (m, b_scc)]),
                l0,
                u1 + m,
            )
            add_constraint!(
                model,
                lin_expr(Tuple{Float64, VarId}[(1.0, lit.vid), (-m, b_scc)]),
                l1 - m,
                u0,
            )
            set_var_bounds!(model, lit.vid, new_lb, new_ub)
            nrewritten += 1
        end

        substitute_scc_by_new_var!(propagator, scc_vid, b_scc)
    end

    return nrewritten
end
