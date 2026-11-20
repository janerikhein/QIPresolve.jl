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

function _quad_expr_state_signature(qe::QuadExpr)
    var_ids = sort!(collect(vars(qe)))
    lin_terms = Tuple{VarId, Float64}[]
    quad_terms = Tuple{VarId, VarId, Float64}[]

    for (index, vid_i) in enumerate(var_ids)
        lin_coeff = get_lin_coeff(qe, vid_i)
        lin_coeff == 0.0 || push!(lin_terms, (vid_i, lin_coeff))

        for vid_j in @view var_ids[index:end]
            quad_coeff = get_quad_coeff(qe, vid_i, vid_j)
            quad_coeff == 0.0 || push!(quad_terms, (vid_i, vid_j, quad_coeff))
        end
    end

    return (
        constant = qe.constant,
        lin_terms = lin_terms,
        quad_terms = quad_terms,
    )
end

function _constraint_state_signature(con::Constraint)
    return (
        id = con.id,
        lhs = con.lhs,
        rhs = con.rhs,
        qe = _quad_expr_state_signature(con.qe),
    )
end

function _model_state_signature(model::QPModel)
    var_terms = sort!(collect((var_id, var.lb, var.ub) for (var_id, var) in model.vars))
    return (
        vars = var_terms,
        cons = [_constraint_state_signature(con) for con in model.cons],
        obj = _quad_expr_state_signature(model.obj_expr),
        infeasible = model.infeasible,
    )
end

"""
    get_builders(model, con)

Build the parity builders induced by an equality constraint.

Construct the mod-2 builder for `con` and, when requested and every non-binary
linear coefficient permits it, also construct the corresponding mod-4
relaxation.

# Arguments
- `model`: Quadratic model that owns the variables referenced by `con`.
- `con`: Equality constraint to translate into parity form.

# Returns
- A tuple `(builder_mod2, builder_mod4)`.
- `builder_mod2` is always an `XorConstraintBuilder`.
- `builder_mod4` is either an `XorConstraintBuilder` or `nothing` when the
  mod-4 relaxation is disabled or not applicable.

# Notes
- The returned builders reflect only parity-relevant coefficients.
- This function assumes `con` is already normalized to integral parity data.

# See also
- [`build_parity_model`](@ref)
"""
function get_builders(model::QPModel, con::Constraint; include_mod4::Bool = true)
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

    (!include_mod4 || discard_mod_4) && return builder_mod2, nothing

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
    build_parity_model(model; parity_strategy=:full)

Build the parity-system view of a quadratic model.

Scan equality constraints in `model`, derive their parity builders, count active
variable participation, and assemble a `ParityModel` ordered by descending
occurrence count. The `parity_strategy` controls whether only mod-2 rows are
generated or whether applicable mod-4 rows are included.

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
function build_parity_model(
    model::QPModel;
    parity_strategy = DEFAULT_PRESOLVE_PARITY_STRATEGY,
)
    parity_strategy = _normalize_parity_strategy(parity_strategy)
    include_mod4 = _parity_includes_mod4(parity_strategy)
    con_builders = XorConstraintBuilder[]

    for con in model.cons
        !is_equality(con) && continue
        builder_mod2, builder_mod4 = get_builders(model, con; include_mod4 = include_mod4)
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

"""
    parity_presolve_phase!(model, propagator, postsolver=nothing; parity_strategy=:full)

Execute one parity presolve phase on a quadratic model.

    Normalize `model` without aggregating parallel constraints, build and
    propagate its parity system, eliminate parity rows, apply any parity-based
    variable fixes or pattern rewrites discovered during the phase, and normalize
    again without aggregation.

# Arguments
- `model`: Quadratic model mutated in place.
- `propagator`: Propagation manager reused across phases.
- `postsolver`: Optional `ParityPostsolver` updated to preserve reconstruction
  data for later postsolve.
- `parity_strategy`: One of `:mod2_basic`, `:mod4_basic`, or `:full`.

# Returns
- Whether the phase changed the model.

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
    postsolver::Union{Nothing, ParityPostsolver} = nothing;
    parity_strategy = DEFAULT_PRESOLVE_PARITY_STRATEGY,
)
    parity_strategy = _normalize_parity_strategy(parity_strategy)
    model.infeasible && return false

    before_normalization = _model_state_signature(model)
    normalize!(model, postsolver; aggregate_parallel = false)
    normalization_changed = before_normalization != _model_state_signature(model)
    model.infeasible && return normalization_changed

    if isempty(model.vars)
        finalize_phase!(propagator)
        return normalization_changed
    end

    parity_model = build_parity_model(model; parity_strategy = parity_strategy)
    if isempty(parity_model.pos_to_var_id) || isempty(parity_model.cons)
        # A row with no parity variables can still assert 0 = 1.
        cleanup!(parity_model)
        model.infeasible = parity_model.infeasible
        finalize_phase!(propagator)
        return normalization_changed
    end

    ensure_literals!(propagator, parity_model.pos_to_var_id)
    propagate!(parity_model, propagator; parity_strategy = parity_strategy)
    if parity_model.infeasible
        model.infeasible = true
        return normalization_changed
    end

    while has_unpivoted_con(parity_model)
        if has_unpivoted_xor_con(parity_model)
            gauss_jordan_xor!(parity_model)
            if parity_model.infeasible
                model.infeasible = true
                return normalization_changed
            end
            propagate!(parity_model, propagator; parity_strategy = parity_strategy)
            if parity_model.infeasible
                model.infeasible = true
                return normalization_changed
            end
            substitute_pivots_in_conjunctive_terms!(parity_model)
            substitute_parity_pivots!(parity_model)
        else
            gauss_jordan_xor_and!(parity_model)
            if parity_model.infeasible
                model.infeasible = true
                return normalization_changed
            end
            propagate!(parity_model, propagator; parity_strategy = parity_strategy)
            if parity_model.infeasible
                model.infeasible = true
                return normalization_changed
            end
        end
    end

    parities_fixed = fix_parities!(model, propagator, postsolver)
    patterns_rewritten = fix_parity_patterns!(model, propagator, postsolver)

    before_post_normalization = _model_state_signature(model)
    normalize!(model, postsolver; aggregate_parallel = false)
    changed = normalization_changed || parities_fixed || patterns_rewritten ||
        before_post_normalization != _model_state_signature(model)

    model.infeasible || finalize_phase!(propagator)
    return changed
end

"""
    parity_presolve!(model, postsolver=nothing; parity_strategy=:full, max_parity_phases=100)
    parity_presolve!(model, propagator, postsolver=nothing; parity_strategy=:full, max_parity_phases=100)

Run parity presolve up to the configured phase limit.

Stop when a phase reports no changes, detects infeasibility, or reaches the phase
limit. A fresh propagation manager is allocated unless `propagator` is supplied.

# Arguments
- `model`: Quadratic model mutated in place.
- `propagator`: Optional propagation manager reused across phases.
- `postsolver`: Optional `ParityPostsolver` that records reconstruction data
  across all phases.
- `max_parity_phases`: Maximum number of phases, including the convergence check. Zero skips parity.
- `parity_strategy`: One of `:mod2_basic`, `:mod4_basic`, or `:full`.

# Returns
- A named tuple with fields `changed`, `domains_changed`,
  `coefficient_changed_constraint_ids`, `phases_run`, and `infeasible`.

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
function parity_presolve!(
    model::QPModel,
    postsolver::Union{Nothing, ParityPostsolver} = nothing;
    parity_strategy = DEFAULT_PRESOLVE_PARITY_STRATEGY,
    max_parity_phases::Integer = DEFAULT_MAX_PARITY_PHASES,
)
    return parity_presolve!(
        model,
        PropagationManager(VarId[]),
        postsolver;
        parity_strategy = parity_strategy,
        max_parity_phases = max_parity_phases,
    )
end

function parity_presolve!(
    model::QPModel,
    propagator::PropagationManager,
    postsolver::Union{Nothing, ParityPostsolver} = nothing;
    parity_strategy = DEFAULT_PRESOLVE_PARITY_STRATEGY,
    max_parity_phases::Integer = DEFAULT_MAX_PARITY_PHASES,
)
    parity_strategy = _normalize_parity_strategy(parity_strategy)
    max_parity_phases >= 0 || throw(ArgumentError("max_parity_phases must be nonnegative"))
    if model.infeasible || iszero(max_parity_phases)
        return (
            changed = false,
            domains_changed = false,
            coefficient_changed_constraint_ids = Set{Int}(),
            infeasible = model.infeasible,
            phases_run = 0,
        )
    end

    before_state = _model_state_signature(model)
    before_domains = _var_domain_snapshot(model)
    before_signatures = _constraint_coefficient_signatures(model)
    normalize!(model, postsolver)
    changed = before_state != _model_state_signature(model)

    phases_run = 0
    while !model.infeasible && phases_run < max_parity_phases
        phases_run += 1
        phase_changed = parity_presolve_phase!(model, propagator, postsolver; parity_strategy = parity_strategy)
        phase_changed || break
        changed = true
    end

    domains_changed = before_domains != _var_domain_snapshot(model)
    coefficient_changed_constraint_ids = _changed_coefficient_constraint_ids(before_signatures, model)
    return (
        changed = changed || before_state != _model_state_signature(model) || domains_changed || !isempty(coefficient_changed_constraint_ids),
        domains_changed = domains_changed,
        coefficient_changed_constraint_ids = coefficient_changed_constraint_ids,
        infeasible = model.infeasible,
        phases_run = phases_run,
    )
end
