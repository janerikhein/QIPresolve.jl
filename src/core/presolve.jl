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

function _symmetrize_formulation!(model::QPModel)
    model.infeasible && return model

    for con in model.cons
        _scale_constraint_by_two!(con)
    end

    return model
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

The pass mutates `model` in place, normalizes and symmetrizes constraints,
alternates parity and residue reductions using the improvement gates expected by
the combined presolver, and returns the mutated model together with postsolve
reconstruction data. Keyword defaults are defined in `QIPresolve.PresolveConfig`.
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

    normalize!(model, postsolver)
    model.infeasible && return PresolveResult(model, postsolver)

    _symmetrize_formulation!(model)
    model.infeasible && return PresolveResult(model, postsolver)

    propagator = PropagationManager(VarId[])

    parity_presolve!(model, propagator, postsolver)
    model.infeasible && return PresolveResult(model, postsolver)

    residue_stats = _residue_presolve_pass!(
        model,
        moduli,
        nothing;
        treewidth_threshold = treewidth_threshold,
    )

    while !model.infeasible && residue_stats.tightened_to_equality
        parity_stats = parity_presolve!(model, propagator, postsolver)
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
