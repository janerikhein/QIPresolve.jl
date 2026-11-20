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
    model.infeasible && return (changed = false, tightened_to_equality = false, infeasible = true)
    normalize!(model, postsolver)
    model.infeasible && return (changed = true, tightened_to_equality = false, infeasible = true)
    isempty(moduli) && return (changed = false, tightened_to_equality = false, infeasible = false)

    changed = false
    tightened_to_equality = false
    for con in model.cons
        candidate_con_ids === nothing || con.id in candidate_con_ids || continue

        result = residue_presolve_constraint!(model, con, moduli; treewidth_threshold = treewidth_threshold)
        changed |= result.changed
        tightened_to_equality |= result.tightened_to_equality
        model.infeasible && break
    end

    return (changed = changed, tightened_to_equality = tightened_to_equality, infeasible = model.infeasible)
end

"""
    presolve!(model; enable_parity=true, enable_residue=true, parity_strategy,
              max_parity_phases, residue_strategy, residue_threshold, treewidth_threshold)

Run the combined parity and residue presolve pipeline.

The pass mutates `model` in place, normalizes and symmetrizes constraints,
alternates parity and residue reductions using the improvement gates expected by
the combined presolver, and returns the mutated model together with postsolve
reconstruction data. Keyword defaults are defined in `QIPresolve.PresolveConfig`.

Set `enable_parity=false` or `enable_residue=false` to skip that reduction
entirely, including its strategy validation. Shared
normalization and postsolve tracking still run, even with both reductions off.
All parity rounds share `max_parity_phases` phases (default 100). A limit of
zero skips parity. Residue presolve can run after parity exhausts its budget.
"""
function presolve!(
        model::QPModel;
        enable_parity::Bool = true,
        enable_residue::Bool = true,
        parity_strategy = DEFAULT_PRESOLVE_PARITY_STRATEGY,
        max_parity_phases::Integer = DEFAULT_MAX_PARITY_PHASES,
        residue_strategy::Symbol = DEFAULT_PRESOLVE_RESIDUE_STRATEGY,
        residue_threshold::Integer = DEFAULT_PRESOLVE_RESIDUE_THRESHOLD,
        treewidth_threshold::Integer = DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD,
    )::PresolveResult
    max_parity_phases >= 0 || throw(ArgumentError("max_parity_phases must be nonnegative"))
    enable_parity && (parity_strategy = _normalize_parity_strategy(parity_strategy))
    postsolver = ParityPostsolver(keys(model.vars))
    moduli = enable_residue ? _generate_residue_moduli(residue_strategy, residue_threshold) : Int[]
    model.infeasible && return PresolveResult(model, postsolver)

    normalize!(model, postsolver)
    model.infeasible && return PresolveResult(model, postsolver)

    propagator = PropagationManager(VarId[])

    remaining_phases = enable_parity ? max_parity_phases : 0
    if remaining_phases > 0
        parity_result = parity_presolve!(
            model, propagator, postsolver;
            parity_strategy = parity_strategy, max_parity_phases = remaining_phases,
        )
        remaining_phases -= parity_result.phases_run
    end
    (model.infeasible || !enable_residue) && return PresolveResult(model, postsolver)

    residue_result = _residue_presolve_pass!(
        model,
        moduli,
        nothing,
        postsolver;
        treewidth_threshold = treewidth_threshold,
    )

    while remaining_phases > 0 && !model.infeasible && residue_result.tightened_to_equality
        parity_result = parity_presolve!(
            model,
            propagator,
            postsolver;
            parity_strategy = parity_strategy,
            max_parity_phases = remaining_phases,
        )
        remaining_phases -= parity_result.phases_run
        model.infeasible && break
        parity_result.domains_changed || break

        residue_result = _residue_presolve_pass!(
            model,
            moduli,
            parity_result.coefficient_changed_constraint_ids,
            postsolver;
            treewidth_threshold = treewidth_threshold,
        )
    end

    return PresolveResult(model, postsolver)
end
