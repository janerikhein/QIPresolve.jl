"""
    PresolveResult

Return value for the combined presolve entry point.

# Fields
- `model`: The mutated presolved model.
- `postsolver`: Reconstruction data for mapping reduced solutions back to the
  original variable space.
"""
struct PresolveResult
    model::QPModel
    postsolver::ParityPostsolver
end

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

function _parity_presolve_fixed_point!(
        model::QPModel,
        propagator::PropagationManager,
        postsolver::ParityPostsolver,
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

    return (
        changed = changed || domains_changed || !isempty(coefficient_changed_constraint_ids),
        domains_changed = domains_changed,
        coefficient_changed_constraint_ids = coefficient_changed_constraint_ids,
        fixed_parities = fixed_parities,
        pattern_rewritten_vars = pattern_rewritten_vars,
        infeasible = model.infeasible,
    )
end

function _residue_presolve_pass!(
        model::QPModel,
        moduli::AbstractVector{<:Integer},
        candidate_con_ids::Union{Nothing, Set{Int}};
        treewidth_threshold::Integer,
    )
    changed = false
    tightened_to_equality = false
    processed = 0
    processed_constraint_ids = Int[]

    (model.infeasible || isempty(moduli)) && return (
        changed = false,
        tightened_to_equality = false,
        processed = 0,
        processed_constraint_ids = processed_constraint_ids,
        infeasible = model.infeasible,
    )

    for con in model.cons
        candidate_con_ids === nothing || con.id in candidate_con_ids || continue

        stats = residue_presolve_constraint!(
            model,
            con,
            moduli;
            treewidth_threshold = treewidth_threshold,
        )

        if stats.processed
            processed += 1
            push!(processed_constraint_ids, con.id)
        end
        changed = changed || stats.changed
        tightened_to_equality = tightened_to_equality || stats.tightened_to_equality

        model.infeasible && break
    end

    return (
        changed = changed,
        tightened_to_equality = tightened_to_equality,
        processed = processed,
        processed_constraint_ids = processed_constraint_ids,
        infeasible = model.infeasible,
    )
end

"""
    presolve!(model; residue_strategy, residue_threshold, treewidth_threshold)

Run the combined parity and residue presolve pipeline.

The pass mutates `model` in place, alternates parity and residue reductions
using the improvement gates expected by the combined presolver, and returns the
mutated model together with postsolve reconstruction data. Keyword defaults are
defined in `QIPresolve.PresolveConfig`.
"""
function presolve!(
        model::QPModel;
        residue_strategy::Symbol = DEFAULT_PRESOLVE_RESIDUE_STRATEGY,
        residue_threshold::Integer = DEFAULT_PRESOLVE_RESIDUE_THRESHOLD,
        treewidth_threshold::Integer = DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD,
    )::PresolveResult
    postsolver = ParityPostsolver(keys(model.vars))
    moduli = _generate_residue_moduli(residue_strategy, residue_threshold)
    model.infeasible && return PresolveResult(model, postsolver)

    propagator = PropagationManager(VarId[])

    _parity_presolve_fixed_point!(model, propagator, postsolver)
    model.infeasible && return PresolveResult(model, postsolver)

    residue_stats = _residue_presolve_pass!(
        model,
        moduli,
        nothing;
        treewidth_threshold = treewidth_threshold,
    )

    while !model.infeasible && residue_stats.tightened_to_equality
        parity_stats = _parity_presolve_fixed_point!(model, propagator, postsolver)
        model.infeasible && break
        parity_stats.domains_changed || break

        residue_stats = _residue_presolve_pass!(
            model,
            moduli,
            parity_stats.coefficient_changed_constraint_ids;
            treewidth_threshold = treewidth_threshold,
        )
    end

    return PresolveResult(model, postsolver)
end
