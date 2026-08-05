const _CoefficientSignature = NamedTuple{
    (:lin_terms, :quad_terms),
    Tuple{
        Vector{Tuple{VarId, Float64}},
        Vector{Tuple{VarId, VarId, Float64}},
    },
}

function _var_domain_snapshot(model::QPModel)
    return Dict{VarId, IntVar}(var_id => var for (var_id, var) in model.vars)
end

function _constraint_coefficient_signature(con::Constraint)::_CoefficientSignature
    var_ids = sort!(collect(vars(con.qe)))
    lin_terms = Tuple{VarId, Float64}[]
    quad_terms = Tuple{VarId, VarId, Float64}[]

    for (index, vid_i) in enumerate(var_ids)
        lin_coeff = get_lin_coeff(con.qe, vid_i)
        lin_coeff == 0.0 || push!(lin_terms, (vid_i, lin_coeff))

        for vid_j in @view var_ids[index:end]
            quad_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            quad_coeff == 0.0 || push!(quad_terms, (vid_i, vid_j, quad_coeff))
        end
    end

    return (lin_terms = lin_terms, quad_terms = quad_terms)
end

function _constraint_coefficient_signatures(model::QPModel)
    signatures = Dict{Int, _CoefficientSignature}()
    sizehint!(signatures, length(model.cons))

    for con in model.cons
        signatures[con.id] = _constraint_coefficient_signature(con)
    end

    return signatures
end

function _changed_coefficient_constraint_ids(
        before::Dict{Int, _CoefficientSignature},
        model::QPModel,
    )
    changed_ids = Set{Int}()

    for con in model.cons
        current = _constraint_coefficient_signature(con)
        if !haskey(before, con.id) || before[con.id] != current
            push!(changed_ids, con.id)
        end
    end

    return changed_ids
end

"""
    get_builders(model, con)

Build the parity builders induced by an equality constraint.

Construct the mod-2 builder for `con` and, when every non-binary linear
coefficient permits it, also construct the corresponding mod-4 relaxation.

# Arguments
- `model`: Quadratic model that owns the variables referenced by `con`.
- `con`: Equality constraint to translate into parity form.

# Returns
- A tuple `(builder_mod2, builder_mod4)`.
- `builder_mod2` is always an `XorConstraintBuilder`.
- `builder_mod4` is either an `XorConstraintBuilder` or `nothing` when the
  mod-4 relaxation is not applicable.

# Notes
- The returned builders reflect only parity-relevant coefficients.
- This function assumes `con` is already normalized to integral parity data.

# See also
- [`build_parity_model`](@ref)
"""
function get_builders(model::QPModel, con::Constraint)
    con_vars = vars(con)

    # Build mod-2 XOR constraint and check if mod-4 relaxation is applicable.
    builder_mod2 = XorConstraintBuilder()
    discard_mod_4 = false

    for vid in con_vars
        is_bin = is_binary(model.vars[vid])
        diag = convert(Int, get_quad_coeff(con.qe, vid, vid))
        lin = convert(Int, get_lin_coeff(con.qe, vid))

        # Mod 4 applies only when every non-binary linear coefficient is even.
        !is_bin && mod(lin, 2) == 1 && (discard_mod_4 = true)

        mod(diag + lin, 2) == 1 && add_par!(builder_mod2, vid)
    end

    con_rhs = convert(Int, con.rhs)
    mod(con_rhs, 2) == 1 && negate!(builder_mod2)

    discard_mod_4 && return builder_mod2, nothing

    builder_mod4 = XorConstraintBuilder()

    for (i, vid_i) in enumerate(con_vars)
        is_bin_i = is_binary(model.vars[vid_i])
        diag_i = convert(Int, get_quad_coeff(con.qe, vid_i, vid_i))
        lin_i = convert(Int, get_lin_coeff(con.qe, vid_i))

        if !is_bin_i
            @assert mod(lin_i, 2) == 0
            mod(lin_i ÷ 2, 2) == 1 && add_par!(builder_mod4, vid_i)
        end

        ci = is_bin_i ? diag_i + lin_i : diag_i
        (mod(ci, 4) >> 1) == 1 && add_par!(builder_mod4, vid_i)

        for j in (i + 1):lastindex(con_vars)
            vid_j = con_vars[j]
            is_bin_j = is_binary(model.vars[vid_j])
            diag_j = convert(Int, get_quad_coeff(con.qe, vid_j, vid_j))
            lin_j = convert(Int, get_lin_coeff(con.qe, vid_j))
            bilin_ij = convert(Int, get_quad_coeff(con.qe, vid_i, vid_j))
            @assert mod(bilin_ij, 2) == 0

            cj = is_bin_j ? diag_j + lin_j : diag_j
            ((mod(bilin_ij ÷ 2, 2) == 1) ⊻ (mod(ci, 2) == 1 && mod(cj, 2) == 1)) &&
                add_conj!(builder_mod4, vid_i, vid_j)
        end
    end

    (mod(con_rhs, 4) >> 1) == 1 && negate!(builder_mod4)

    return builder_mod2, builder_mod4
end

"""
    count_builder_occurrences!(var_counts, builder)

Accumulate active variable occurrences from a parity builder.

Count active parity literals and active conjunctive terms in `builder`, adding
their contributions into `var_counts`.

# Arguments
- `var_counts`: Dictionary updated in place with occurrence counts per variable.
- `builder`: Builder whose active entries should be counted.

# Returns
- The mutated `var_counts` dictionary.

# Side Effects
- Mutates `var_counts`.

# Notes
- Disabled builder entries are ignored.

# See also
- [`build_parity_model`](@ref)
"""
function count_builder_occurrences!(
    var_counts::Dict{VarId, Int},
    builder::XorConstraintBuilder,
)
    for (vid, is_active) in builder.par
        is_active || continue
        var_counts[vid] = get(var_counts, vid, 0) + 1
    end

    for ((vid_i, vid_j), is_active) in builder.conj
        is_active || continue
        var_counts[vid_i] = get(var_counts, vid_i, 0) + 1
        var_counts[vid_j] = get(var_counts, vid_j, 0) + 1
    end

    return var_counts
end

"""
    build_parity_model(model)

Build the parity-system view of a quadratic model.

Scan equality constraints in `model`, derive their parity builders, count active
variable participation, and assemble a `ParityModel` ordered by descending
occurrence count.

# Arguments
- `model`: Quadratic model to translate into parity constraints.

# Returns
- A `ParityModel` containing the parity-variable ordering and built constraints.

# Notes
- Only equality constraints contribute rows.
- Variables with no active parity participation are omitted from the returned
  parity model.

# See also
- [`get_builders`](@ref)
- [`count_builder_occurrences!`](@ref)
"""
function build_parity_model(model::QPModel)
    con_builders = XorConstraintBuilder[]

    for con in model.cons
        !is_equality(con) && continue
        builder_mod2, builder_mod4 = get_builders(model, con)
        push!(con_builders, builder_mod2)
        builder_mod4 !== nothing && push!(con_builders, builder_mod4)
    end

    var_counts = Dict{VarId, Int}()
    for builder in con_builders
        count_builder_occurrences!(var_counts, builder)
    end

    pos_to_var_id = collect(keys(var_counts))
    sort!(pos_to_var_id; by = vid -> (-var_counts[vid], vid))
    var_id_to_pos = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var_id))

    nparity_vars = length(pos_to_var_id)
    cons = [
        build(builder, nparity_vars, var_id_to_pos) for builder in con_builders
    ]

    return ParityModel(var_id_to_pos, pos_to_var_id, cons)
end

parity_presolve_phase!(model::QPModel, propagator::PropagationManager) = parity_presolve_phase!(model, propagator, nothing)

"""
    parity_presolve_phase!(model, propagator, postsolver=nothing)

Execute one parity presolve phase on a quadratic model.

Normalize `model`, build and propagate its parity system, eliminate parity rows,
and apply any parity-based variable fixes or pattern rewrites discovered during
the phase.

# Arguments
- `model`: Quadratic model mutated in place.
- `propagator`: Propagation manager reused across phases.
- `postsolver`: Optional `ParityPostsolver` updated to preserve reconstruction
  data for later postsolve.

# Returns
- A named tuple with fields `changed`, `fixed_parities`, and
  `pattern_rewritten_vars`.

# Side Effects
- Mutates `model` and `propagator`.
- May mutate `postsolver` when it is provided.
- May set `model.infeasible = true`.

# Notes
- The phase finalizes the propagator before returning from non-infeasible early
  exits.
- The caller is responsible for repeating phases until a fixed point is reached.

# See also
- [`parity_presolve!`](@ref)
- [`build_parity_model`](@ref)
"""
function parity_presolve_phase!(
    model::QPModel,
    propagator::PropagationManager,
    postsolver::Union{Nothing, ParityPostsolver},
)
    model.infeasible && return (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)

    normalize!(model, postsolver)
    model.infeasible && return (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)

    if isempty(model.vars)
        finalize_phase!(propagator)
        return (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)
    end

    parity_model = build_parity_model(model)
    if isempty(parity_model.pos_to_var_id) || isempty(parity_model.cons)
        finalize_phase!(propagator)
        return (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)
    end

    ensure_literals!(propagator, parity_model.pos_to_var_id)
    propagate!(parity_model, propagator)
    if parity_model.infeasible
        model.infeasible = true
        return (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)
    end

    while has_unpivoted_con(parity_model)
        if has_unpivoted_xor_con(parity_model)
            gauss_jordan_xor!(parity_model)
            propagate!(parity_model, propagator)
            if parity_model.infeasible
                model.infeasible = true
                return (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)
            end
            substitute_pivots_in_conjunctive_terms!(parity_model)
            substitute_parity_pivots!(parity_model)
        else
            gauss_jordan_xor_and!(parity_model)
            propagate!(parity_model, propagator)
            if parity_model.infeasible
                model.infeasible = true
                return (changed = false, fixed_parities = 0, pattern_rewritten_vars = 0)
            end
        end
    end

    parities_fixed = fix_parities!(model, propagator, postsolver)
    pattern_rewritten_vars = fix_parity_patterns!(model, propagator, postsolver)
    changed = parities_fixed > 0 || pattern_rewritten_vars > 0

    finalize_phase!(propagator)
    return (changed = changed, fixed_parities = parities_fixed, pattern_rewritten_vars = pattern_rewritten_vars)
end

parity_presolve!(model::QPModel) = parity_presolve!(model, nothing)

parity_presolve!(
    model::QPModel,
    propagator::PropagationManager,
) = parity_presolve!(model, propagator, nothing)

"""
    parity_presolve!(model, postsolver=nothing)
    parity_presolve!(model, propagator, postsolver=nothing)

Run parity presolve to a fixed point.

Repeatedly execute parity presolve phases until `model` becomes infeasible or a
phase reports no further parity-driven changes. The two-argument form allocates a
fresh propagation manager; the three-argument form reuses `propagator`.

# Arguments
- `model`: Quadratic model mutated in place.
- `propagator`: Optional propagation manager reused across phases.
- `postsolver`: Optional `ParityPostsolver` that records reconstruction data
  across all phases.

# Returns
- A named tuple with fields `changed`, `domains_changed`,
  `coefficient_changed_constraint_ids`, `fixed_parities`,
  `pattern_rewritten_vars`, and `infeasible`.

# Side Effects
- Mutates `model`.
- Allocates and mutates an internal `PropagationManager` unless one is supplied.
- May mutate `postsolver` when it is provided.

# Notes
- The returned `changed` field also reports normalization-only domain or
  coefficient changes observed across the full fixed-point run.

# See also
- [`parity_presolve_phase!`](@ref)
"""
function parity_presolve!(model::QPModel, postsolver::Union{Nothing, ParityPostsolver})
    return parity_presolve!(model, PropagationManager(VarId[]), postsolver)
end

function parity_presolve!(
    model::QPModel,
    propagator::PropagationManager,
    postsolver::Union{Nothing, ParityPostsolver},
)
    before_domains = _var_domain_snapshot(model)
    before_signatures = _constraint_coefficient_signatures(model)

    changed = false
    fixed_parities = 0
    pattern_rewritten_vars = 0

    while true
        stats = parity_presolve_phase!(model, propagator, postsolver)
        changed = changed || stats.changed
        fixed_parities += stats.fixed_parities
        pattern_rewritten_vars += stats.pattern_rewritten_vars

        (model.infeasible || !stats.changed) && break
    end

    domains_changed = before_domains != _var_domain_snapshot(model)
    coefficient_changed_constraint_ids =
        _changed_coefficient_constraint_ids(before_signatures, model)

    #add_binary_implications!(model, propagator)
    return (
        changed = changed || domains_changed || !isempty(coefficient_changed_constraint_ids),
        domains_changed = domains_changed,
        coefficient_changed_constraint_ids = coefficient_changed_constraint_ids,
        fixed_parities = fixed_parities,
        pattern_rewritten_vars = pattern_rewritten_vars,
        infeasible = model.infeasible,
    )
end
