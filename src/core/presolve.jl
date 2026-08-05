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

function _residue_presolve_pass!(
        model::QPModel,
        moduli::AbstractVector{<:Integer},
        candidate_con_ids::Union{Nothing, Set{Int}},
        postsolver::Union{Nothing, ParityPostsolver} = nothing;
        treewidth_threshold::Integer,
    )
    changed = false
    tightened_to_equality = false
    processed = 0
    processed_constraint_ids = Int[]

    model.infeasible && return (
        changed = false,
        tightened_to_equality = false,
        processed = 0,
        processed_constraint_ids = processed_constraint_ids,
        infeasible = model.infeasible,
    )

    normalize!(model, postsolver)
    model.infeasible && return (
        changed = true,
        tightened_to_equality = false,
        processed = 0,
        processed_constraint_ids = processed_constraint_ids,
        infeasible = true,
    )

    isempty(moduli) && return (
        changed = false,
        tightened_to_equality = false,
        processed = 0,
        processed_constraint_ids = processed_constraint_ids,
        infeasible = false,
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

    propagator = PropagationManager(VarId[])

    parity_presolve!(model, propagator, postsolver)
    model.infeasible && return PresolveResult(model, postsolver)

    residue_stats = _residue_presolve_pass!(
        model,
        moduli,
        nothing,
        postsolver;
        treewidth_threshold = treewidth_threshold,
    )

    while !model.infeasible && residue_stats.tightened_to_equality
        parity_stats = parity_presolve!(model, propagator, postsolver)
        model.infeasible && break
        parity_stats.domains_changed || break

        residue_stats = _residue_presolve_pass!(
            model,
            moduli,
            parity_stats.coefficient_changed_constraint_ids,
            postsolver;
            treewidth_threshold = treewidth_threshold,
        )
    end

    return PresolveResult(model, postsolver)
end
