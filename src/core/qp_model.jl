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
- Whether any variable was transformed.

# Side Effects
- Mutates `model`.
- May mutate `postsolver`.
"""
fix_parities!(model::QPModel, propagator::PropagationManager) = fix_parities!(model, propagator, nothing)

# Each participating variable is replaced once by 2q+b or 2q+1-b. Its
# expansion has absolute coefficient sum at most 4 (16 for a product).
# Bound the intermediates before changing any member of an unbounded group.
function _parity_rewrite_fits(qe::QuadExpr, ids, finite_bounds...)
    any(id -> has_var(qe, id), ids) || return true
    limit = maxintfloat() / 16
    largest = 0.0
    count = 0
    for value in Iterators.flatten((finite_bounds, (qe.constant,), lin(qe), quad(qe)))
        iszero(value) && continue
        isfinite(value) || return false
        count += 1
        largest = max(largest, abs(value))
        if !isinteger(value)
            # Float64 coefficients are dyadic. Allow for their fractional bits
            # when bounding exact additions, including shifts of a row bound.
            den = denominator(Rational{BigInt}(value))
            limit = min(limit, ldexp(maxintfloat() / 16, -trailing_zeros(den)))
        end
    end
    return count == 0 || largest < limit / count
end

function _unbounded_parity_rewrite_fits(model::QPModel, ids)
    for id in ids
        var = model.vars[id]
        var.lb <= var.ub && var.lb < Inf && var.ub > -Inf || return false
        # Leave room for rounding each finite endpoint and forming linking rows.
        for endpoint in (var.lb, var.ub)
            isfinite(endpoint) && abs(endpoint) >= maxintfloat() / 4 && return false
        end
    end
    _parity_rewrite_fits(model.obj_expr, ids) || return false
    for con in model.cons
        lower = isfinite(con.lhs) ? abs(con.lhs) : 0.0
        upper = isfinite(con.rhs) ? abs(con.rhs) : 0.0
        _parity_rewrite_fits(con.qe, ids, lower, upper) || return false
    end
    return true
end

function fix_parities!(
    model::QPModel,
    propagator::PropagationManager,
    postsolver::Union{Nothing, ParityPostsolver},
)
    changed = false

    for lit in propagator.pos_to_lit
        lit.neg && continue
        haskey(model.vars, lit.vid) || continue

        val = fixed_value(propagator, lit.vid)
        val === nothing && continue

        offset = val ? 1.0 : 0.0
        var = model.vars[lit.vid]
        if !isfinite(var.lb) || !isfinite(var.ub)
            _unbounded_parity_rewrite_fits(model, (lit.vid,)) || continue
        end
        new_lb = ceil((ceil(var.lb) - offset) / 2.0)
        new_ub = floor((floor(var.ub) - offset) / 2.0)
        new_lb == 0.0 && (new_lb = 0.0)
        new_ub == 0.0 && (new_ub = 0.0)

        postsolver !== nothing && append_fixed_bit!(postsolver, lit.vid, val)
        affine_transform!(model, lit.vid, 2.0, offset)
        set_var_bounds!(model, lit.vid, new_lb, new_ub)
        changed = true
    end

    return changed
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
- Whether any variable was rewritten.

# Side Effects
- Mutates `model` and `propagator`.
- May mutate `postsolver`.
"""
function fix_parity_patterns!(
    model::QPModel,
    propagator::ParityPropagator,
    postsolver::Union{Nothing, ParityPostsolver},
)
    propagator.nreps == length(propagator.parent_pos) && return false
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

    changed = false

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
        if any(lit -> !isfinite(model.vars[lit.vid].lb) || !isfinite(model.vars[lit.vid].ub), component_lits)
            _unbounded_parity_rewrite_fits(model, [lit.vid for lit in component_lits]) || continue
        end
        @assert length(unique(lit.vid for lit in component_lits)) == length(component_lits)

        b_scc = add_var!(model, IntVar(0.0, 1.0))
        postsolver !== nothing && ensure_tracked_var!(postsolver, b_scc)

        for lit in component_lits
            var = model.vars[lit.vid]
            e_lo = ceil(ceil(var.lb) / 2.0)
            e_hi = floor(floor(var.ub) / 2.0)
            o_lo = ceil((ceil(var.lb) - 1.0) / 2.0)
            o_hi = floor((floor(var.ub) - 1.0) / 2.0)

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

            new_lb = min(l0, l1)
            new_ub = max(u0, u1)
            new_lb == 0.0 && (new_lb = 0.0)
            new_ub == 0.0 && (new_ub = 0.0)

            if isfinite(new_lb) && isfinite(new_ub)
                m = max(0.0, u0 - l1, u1 - l0)
                add_constraint!(model, lin_expr([(1.0, lit.vid), (m, b_scc)]), l0, u1 + m)
                add_constraint!(model, lin_expr([(1.0, lit.vid), (-m, b_scc)]), l1 - m, u0)
            elseif l0 != l1
                add_constraint!(model, lin_expr([(1.0, lit.vid), (l0 - l1, b_scc)]), l0, Inf)
            elseif u0 != u1
                add_constraint!(model, lin_expr([(1.0, lit.vid), (u0 - u1, b_scc)]), -Inf, u0)
            end
            set_var_bounds!(model, lit.vid, new_lb, new_ub)
        end

        substitute_scc_by_new_var!(propagator, scc_vid, b_scc)
        changed = true
    end

    return changed
end
