@inline _canonicalize_zero(x::Float64) = (x == 0.0 ? 0.0 : x)

@inline function _is_shifted_binary(var::IntVar)
    return isinteger(var.lb) && isinteger(var.ub) && var.ub == var.lb + 1.0 && var.lb != 0.0
end

"""
    normalize!(qe)

Remove inactive variables from `qe`.

Delete variables whose linear and quadratic coefficients are all zero and refresh
the cached integrality flag.
"""
function normalize!(qe::QuadExpr)
    for pos in reverse(1:qe.nvars)
        phys_pos = qe.perm[pos]

        qe.lin_buf[phys_pos] != 0 && continue

        all_zero = true
        for other_pos in 1:qe.nvars
            phys_other = qe.perm[other_pos]
            q = if phys_pos <= phys_other
                @inbounds qe.quad_buf[phys_pos, phys_other]
            else
                @inbounds qe.quad_buf[phys_other, phys_pos]
            end
            if q != 0
                all_zero = false
                break
            end
        end
        all_zero && remove_var!(qe, qe.pos_to_var[pos])
    end

    return _refresh_is_integer!(qe)
end

"""
    normalize!(con::Constraint) -> Constraint

Move the constant term from the quadratic expression into the bounds.

After normalization, `con.qe.constant == 0.0` and `lhs`/`rhs` are shifted by the
previous constant so the constraint is equivalent. Integral expressions also
tighten finite bounds to integers.
"""
function normalize!(con::Constraint)
    normalize!(con.qe)
    con.lhs -= con.qe.constant
    con.rhs -= con.qe.constant
    con.qe.constant = 0.0

    if isinteger(con.qe)
        isfinite(con.lhs) && (con.lhs = ceil(con.lhs))
        isfinite(con.rhs) && (con.rhs = floor(con.rhs))
    end

    con.lhs = _canonicalize_zero(con.lhs)
    con.rhs = _canonicalize_zero(con.rhs)
    return con
end

"""
    symmetrize!(con::Constraint) -> Constraint

Symmetrize the quadratic form in-place.

For the dense-buffer representation, this replaces `Q` by `Q + Q'` and rescales
the linear term and bounds so the constraint remains equivalent:

    lhs <= x'Qx + c'x + constant <= rhs

becomes

    2lhs <= x'(Q+Q')x + 2c'x + 2constant <= 2rhs
"""
function symmetrize!(con::Constraint)
    quad_mat = quad(con.qe)
    lin_vec = lin(con.qe)

    quad_mat .+= transpose(quad_mat)
    lin_vec .*= 2
    con.qe.constant *= 2
    con.lhs *= 2
    con.rhs *= 2
    _refresh_is_integer!(con.qe)
    return con
end

"""
    scale_gcd!(con::Constraint) -> Bool

Normalize `con`, divide all coefficients by their greatest common divisor, and
tighten finite bounds to integers. Returns `true` iff the constraint was scaled.
Constraints with non-integral data are left unchanged after normalization.
"""
function scale_gcd!(con::Constraint)
    normalize!(con)
    is_integer(con) || return false

    var_ids = collect(vars(con.qe))
    nvars = length(var_ids)
    g = 0

    for i in 1:nvars
        vid_i = var_ids[i]

        lin_coeff = get_lin_coeff(con.qe, vid_i)
        if lin_coeff != 0.0
            g = gcd(g, abs(trunc(Int, lin_coeff)))
        end
        g == 1 && return false

        diag_coeff = get_quad_coeff(con.qe, vid_i, vid_i)
        if diag_coeff != 0.0
            g = gcd(g, abs(trunc(Int, diag_coeff)))
        end
        g == 1 && return false

        for j in (i + 1):nvars
            vid_j = var_ids[j]
            bilinear_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            bilinear_coeff == 0.0 && continue
            g = gcd(g, abs(trunc(Int, bilinear_coeff) ÷ 2))
            g == 1 && return false
        end
    end

    g <= 1 && return false

    scale = Float64(g)

    for vid in var_ids
        lin_coeff = get_lin_coeff(con.qe, vid)
        lin_coeff == 0.0 || set_lin_coeff!(con.qe, vid, lin_coeff / scale)

        diag_coeff = get_quad_coeff(con.qe, vid, vid)
        diag_coeff == 0.0 || set_quad_coeff!(con.qe, vid, vid, diag_coeff / scale)
    end

    for i in 1:nvars
        vid_i = var_ids[i]
        for j in (i + 1):nvars
            vid_j = var_ids[j]
            bilinear_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            bilinear_coeff == 0.0 || set_quad_coeff!(con.qe, vid_i, vid_j, bilinear_coeff / scale)
        end
    end

    if isfinite(con.lhs)
        con.lhs = _canonicalize_zero(ceil(con.lhs / scale))
    end
    if isfinite(con.rhs)
        con.rhs = _canonicalize_zero(floor(con.rhs / scale))
    end

    return true
end

"""
    scale_constraints_gcd!(model::QPModel) -> Int

Scale all integral constraints in `model` by the gcd of their coefficients and
tighten finite bounds accordingly. Returns the number of constraints scaled and
marks the model infeasible if a scaled constraint gets inconsistent bounds.
"""
function scale_constraints_gcd!(model::QPModel)
    model.infeasible && return 0

    nscaled = 0
    for con in model.cons
        scale_gcd!(con) && (nscaled += 1)
        if con.lhs > con.rhs
            model.infeasible = true
            return nscaled
        end
    end

    return nscaled
end

function _tighten_singleton_constraint!(model::QPModel, con::Constraint)
    var_id = vars(con.qe)[1]
    coeff = get_lin_coeff(con.qe, var_id)
    @assert coeff != 0.0

    var_bounds = model.vars[var_id]
    bound1 = con.lhs / coeff
    bound2 = con.rhs / coeff
    singleton_lb = ceil(min(bound1, bound2))
    singleton_ub = floor(max(bound1, bound2))
    new_lb = max(var_bounds.lb, singleton_lb)
    new_ub = min(var_bounds.ub, singleton_ub)
    new_lb = _canonicalize_zero(new_lb)
    new_ub = _canonicalize_zero(new_ub)

    if new_lb > new_ub
        model.infeasible = true
        return false
    end

    if new_lb == var_bounds.lb && new_ub == var_bounds.ub
        return false
    end

    model.vars[var_id] = IntVar(new_lb, new_ub)
    return true
end

function _fix_vars_and_singletons!(model::QPModel, postsolver::Union{Nothing, ParityPostsolver})
    model.infeasible && return false

    changed = false
    fixed_var_ids = VarId[]

    for (var_id, var) in model.vars
        if var.lb > var.ub
            model.infeasible = true
            return changed
        end
        var.lb == var.ub || continue
        push!(fixed_var_ids, var_id)
    end

    for var_id in fixed_var_ids
        haskey(model.vars, var_id) || continue
        var = model.vars[var_id]
        postsolver !== nothing && register_fixed_var!(postsolver, var_id, var.lb)
        var_bound_shift!(model, var_id, var.lb)
        @assert model.vars[var_id].lb == model.vars[var_id].ub == 0.0
        for con in model.cons
            remove_var!(con, var_id)
        end
        remove_var!(model.obj_expr, var_id)
        delete!(model.vars, var_id)
        changed = true
    end

    normalize!(model.obj_expr)

    for i in reverse(eachindex(model.cons))
        con = model.cons[i]
        normalize!(con)

        if is_empty(con.qe)
            deleteat!(model.cons, i)
            changed = true
            continue
        end

        if is_singleton(con.qe)
            tightened = _tighten_singleton_constraint!(model, con)
            model.infeasible && return true
            changed = changed || tightened
            deleteat!(model.cons, i)
            changed = true
        end
    end

    return changed
end

function _normalize_shifted_binaries!(model::QPModel, postsolver::Union{Nothing, ParityPostsolver})
    model.infeasible && return false

    changed = false
    shifted_binary_ids = VarId[]

    for (var_id, var) in model.vars
        _is_shifted_binary(var) || continue
        push!(shifted_binary_ids, var_id)
    end

    for var_id in shifted_binary_ids
        haskey(model.vars, var_id) || continue
        var = model.vars[var_id]
        _is_shifted_binary(var) || continue

        postsolver !== nothing && add_reconstruction_offset!(postsolver, var_id, var.lb)
        affine_transform!(model, var_id, 1.0, var.lb)
        set_var_bounds!(model, var_id, 0.0, 1.0)
        changed = true
    end

    return changed
end

function _fold_binary_diagonal!(qe::QuadExpr, binary_var_ids::AbstractVector{VarId})
    changed = false

    for var_id in binary_var_ids
        coeff = get_quad_coeff(qe, var_id, var_id)
        coeff == 0.0 && continue
        add_lin_coeff!(qe, var_id, coeff) || continue
        set_quad_coeff!(qe, var_id, var_id, 0.0)
        changed = true
    end

    changed && normalize!(qe)
    return changed
end

function _fold_binary_diagonal!(con::Constraint, binary_var_ids::AbstractVector{VarId})
    changed = _fold_binary_diagonal!(con.qe, binary_var_ids)
    changed || return false
    normalize!(con)
    return true
end

function _fold_binary_diagonal!(model::QPModel)
    isempty(model.vars) && return false

    binary_var_ids = [var_id for (var_id, var) in model.vars if is_binary(var)]
    isempty(binary_var_ids) && return false

    changed = _fold_binary_diagonal!(model.obj_expr, binary_var_ids)

    for con in model.cons
        changed = _fold_binary_diagonal!(con, binary_var_ids) || changed
    end

    return changed
end

@inline _is_integer_after_doubling(coeff::Float64) = isinteger(2.0 * coeff)

function _is_even_integer_after_doubling(coeff::Float64)
    doubled = 2.0 * coeff
    return isinteger(doubled) && iseven(trunc(Int, doubled))
end

function _can_scale_by_two_to_integer(con::Constraint)
    var_ids = collect(vars(con.qe))

    for (idx, vid_i) in enumerate(var_ids)
        _is_integer_after_doubling(get_lin_coeff(con.qe, vid_i)) || return false
        _is_integer_after_doubling(get_quad_coeff(con.qe, vid_i, vid_i)) || return false

        for vid_j in @view var_ids[(idx + 1):end]
            _is_even_integer_after_doubling(get_quad_coeff(con.qe, vid_i, vid_j)) ||
                return false
        end
    end

    return true
end

function _scale_constraint_by_two!(con::Constraint)
    var_ids = collect(vars(con.qe))

    for vid in var_ids
        lin_coeff = get_lin_coeff(con.qe, vid)
        lin_coeff == 0.0 || set_lin_coeff!(con.qe, vid, 2.0 * lin_coeff)

        diag_coeff = get_quad_coeff(con.qe, vid, vid)
        diag_coeff == 0.0 || set_quad_coeff!(con.qe, vid, vid, 2.0 * diag_coeff)
    end

    for (idx, vid_i) in enumerate(var_ids)
        for vid_j in @view var_ids[(idx + 1):end]
            quad_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            quad_coeff == 0.0 || set_quad_coeff!(con.qe, vid_i, vid_j, 2.0 * quad_coeff)
        end
    end

    con.qe.constant *= 2.0
    con.lhs *= 2.0
    con.rhs *= 2.0
    return normalize!(con)
end

function _symmetrize_quadratic_form!(model::QPModel)
    model.infeasible && return 0

    nsymmetrized = 0
    for con in model.cons
        normalize!(con)
        is_integer(con) && continue
        _can_scale_by_two_to_integer(con) || continue

        _scale_constraint_by_two!(con)
        nsymmetrized += 1
        if con.lhs > con.rhs
            model.infeasible = true
            return nsymmetrized
        end
    end

    return nsymmetrized
end

function _constraint_coefficient_key(con::Constraint)
    var_ids = sort!(collect(vars(con.qe)))
    lin_terms = Tuple{VarId, Float64}[]
    quad_terms = Tuple{VarId, VarId, Float64}[]

    for (idx, vid_i) in enumerate(var_ids)
        lin_coeff = get_lin_coeff(con.qe, vid_i)
        lin_coeff == 0.0 || push!(lin_terms, (vid_i, lin_coeff))

        for vid_j in @view var_ids[idx:end]
            quad_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            quad_coeff == 0.0 || push!(quad_terms, (vid_i, vid_j, quad_coeff))
        end
    end

    return (lin_terms = lin_terms, quad_terms = quad_terms)
end

function _negated_constraint_key(key)
    lin_terms = [(var_id, -coeff) for (var_id, coeff) in key.lin_terms]
    quad_terms = [(var_i, var_j, -coeff) for (var_i, var_j, coeff) in key.quad_terms]
    return (lin_terms = lin_terms, quad_terms = quad_terms)
end

function _tighten_representative_bounds!(rep::Constraint, lhs::Float64, rhs::Float64)
    old_lhs = rep.lhs
    old_rhs = rep.rhs
    rep.lhs = _canonicalize_zero(max(rep.lhs, lhs))
    rep.rhs = _canonicalize_zero(min(rep.rhs, rhs))
    return rep.lhs != old_lhs || rep.rhs != old_rhs
end

function aggregate_parallel_constraints!(model::QPModel)
    model.infeasible && return 0

    representatives = Constraint[]
    signature_to_index = Dict{Any, Int}()
    changed = false

    for con in model.cons
        normalize!(con)
        key = _constraint_coefficient_key(con)
        negated_key = _negated_constraint_key(key)

        if haskey(signature_to_index, key)
            rep = representatives[signature_to_index[key]]
            changed = _tighten_representative_bounds!(rep, con.lhs, con.rhs) || changed
            changed = true
        elseif haskey(signature_to_index, negated_key)
            rep = representatives[signature_to_index[negated_key]]
            changed = _tighten_representative_bounds!(rep, -con.rhs, -con.lhs) || changed
            changed = true
        else
            signature_to_index[key] = length(representatives) + 1
            push!(representatives, con)
        end

        if !isempty(representatives) && last(representatives).lhs > last(representatives).rhs
            model.infeasible = true
            break
        end
        if haskey(signature_to_index, key)
            rep = representatives[signature_to_index[key]]
            if rep.lhs > rep.rhs
                model.infeasible = true
                break
            end
        elseif haskey(signature_to_index, negated_key)
            rep = representatives[signature_to_index[negated_key]]
            if rep.lhs > rep.rhs
                model.infeasible = true
                break
            end
        end
    end

    if changed || length(representatives) != length(model.cons)
        model.cons = representatives
        return 1
    end

    return 0
end

normalize!(model::QPModel; kwargs...) = normalize!(model, nothing; kwargs...)

"""
    normalize!(model, postsolver=nothing; kwargs...)

Normalize `model` to a presolve-stable form.

Keyword arguments control the normalization steps. By default all steps run:
fixed-variable and singleton cleanup, shifted-binary standardization, binary
diagonal folding, idempotent quadratic-form symmetrization, GCD scaling, and
parallel-constraint aggregation.
"""
function normalize!(
        model::QPModel,
        postsolver::Union{Nothing, ParityPostsolver};
        remove_fixed_vars::Bool = true,
        standardize_binaries::Bool = true,
        fold_binary_diagonal::Bool = true,
        symmetrize_quadratic::Bool = true,
        scale_gcd::Bool = true,
        aggregate_parallel::Bool = true,
    )
    model.infeasible && return model

    while true
        changed = false

        if remove_fixed_vars
            changed = _fix_vars_and_singletons!(model, postsolver) || changed
            model.infeasible && return model
        end

        if standardize_binaries
            changed = _normalize_shifted_binaries!(model, postsolver) || changed
            model.infeasible && return model
        end

        if fold_binary_diagonal
            changed = _fold_binary_diagonal!(model) || changed
            model.infeasible && return model
        end

        if symmetrize_quadratic
            changed = (_symmetrize_quadratic_form!(model) > 0) || changed
            model.infeasible && return model
        end

        if scale_gcd
            changed = (scale_constraints_gcd!(model) > 0) || changed
            model.infeasible && return model
        end

        if aggregate_parallel
            changed = (aggregate_parallel_constraints!(model) > 0) || changed
            model.infeasible && return model
        end

        changed || break
    end

    return model
end

fix_vars!(model::QPModel) = fix_vars!(model, nothing)

"""
    fix_vars!(model, postsolver=nothing)

Apply one fixed-variable and singleton cleanup pass to `model`.
"""
function fix_vars!(model::QPModel, postsolver::Union{Nothing, ParityPostsolver})
    _fix_vars_and_singletons!(model, postsolver)
    return model
end
