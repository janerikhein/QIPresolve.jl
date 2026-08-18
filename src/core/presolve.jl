"""
    PresolveResult

Return value for the combined presolve entry point.

# Fields
- `model`: The mutated presolved model.
- `postsolver`: Reconstruction data for mapping reduced solutions back to the
  original variable space.
- `parity_stats`: Aggregate statistics from parity presolve.
- `residue_stats`: Aggregate statistics from residue presolve.
"""
struct PresolveResult
    model::QPModel
    postsolver::ParityPostsolver
    parity_stats::ParityStats
    residue_stats::ResidueStats
end

PresolveResult(model::QPModel, postsolver::ParityPostsolver) =
    PresolveResult(model, postsolver, ParityStats(), ResidueStats())

PresolveResult(model::QPModel, postsolver::ParityPostsolver, parity_stats::ParityStats) =
    PresolveResult(model, postsolver, parity_stats, ResidueStats())

function _residue_presolve_pass!(
        model::QPModel,
        moduli::AbstractVector{<:Integer},
        candidate_con_ids::Union{Nothing, Set{Int}},
        postsolver::Union{Nothing, ParityPostsolver} = nothing;
        treewidth_threshold::Integer,
    )
    return _residue_presolve_pass!(
        model,
        moduli,
        candidate_con_ids,
        postsolver,
        _ResidueStatsAccumulator();
        treewidth_threshold = treewidth_threshold,
    )
end

function _residue_presolve_pass!(
        model::QPModel,
        moduli::AbstractVector{<:Integer},
        candidate_con_ids::Union{Nothing, Set{Int}},
        stats_accumulator::_ResidueStatsAccumulator;
        treewidth_threshold::Integer,
    )
    return _residue_presolve_pass!(
        model,
        moduli,
        candidate_con_ids,
        nothing,
        stats_accumulator;
        treewidth_threshold = treewidth_threshold,
    )
end

function _residue_presolve_pass!(
        model::QPModel,
        moduli::AbstractVector{<:Integer},
        candidate_con_ids::Union{Nothing, Set{Int}},
        postsolver::Union{Nothing, ParityPostsolver},
        stats_accumulator::_ResidueStatsAccumulator;
        treewidth_threshold::Integer,
    )
    start_time = time()
    try
        return _residue_presolve_pass_impl!(
            model,
            moduli,
            candidate_con_ids,
            postsolver,
            stats_accumulator;
            treewidth_threshold = treewidth_threshold,
        )
    finally
        stats_accumulator.stats.residue_presolve_time += time() - start_time
    end
end

function _residue_presolve_pass_impl!(
        model::QPModel,
        moduli::AbstractVector{<:Integer},
        candidate_con_ids::Union{Nothing, Set{Int}},
        postsolver::Union{Nothing, ParityPostsolver},
        stats_accumulator::_ResidueStatsAccumulator;
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
        residue_stats = stats_accumulator.stats,
    )

    normalize!(model, postsolver)
    model.infeasible && return (
        changed = true,
        tightened_to_equality = false,
        processed = 0,
        processed_constraint_ids = processed_constraint_ids,
        infeasible = true,
        residue_stats = stats_accumulator.stats,
    )

    isempty(moduli) && return (
        changed = false,
        tightened_to_equality = false,
        processed = 0,
        processed_constraint_ids = processed_constraint_ids,
        infeasible = false,
        residue_stats = stats_accumulator.stats,
    )

    for con in model.cons
        candidate_con_ids === nothing || con.id in candidate_con_ids || continue

        stats = _residue_presolve_constraint_impl!(
            model,
            con,
            moduli,
            stats_accumulator;
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
        residue_stats = stats_accumulator.stats,
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
    parity_stats_accumulator = _ParityStatsAccumulator()
    residue_stats_accumulator = _ResidueStatsAccumulator()
    moduli = _generate_residue_moduli(residue_strategy, residue_threshold)
    model.infeasible && return PresolveResult(
        model,
        postsolver,
        parity_stats_accumulator.stats,
        residue_stats_accumulator.stats,
    )

    normalize!(model, postsolver)
    model.infeasible && return PresolveResult(
        model,
        postsolver,
        parity_stats_accumulator.stats,
        residue_stats_accumulator.stats,
    )

    propagator = PropagationManager(VarId[])

    parity_presolve!(model, propagator, postsolver, parity_stats_accumulator)
    model.infeasible && return PresolveResult(
        model,
        postsolver,
        parity_stats_accumulator.stats,
        residue_stats_accumulator.stats,
    )

    residue_stats = _residue_presolve_pass!(
        model,
        moduli,
        nothing,
        postsolver,
        residue_stats_accumulator;
        treewidth_threshold = treewidth_threshold,
    )

    while !model.infeasible && residue_stats.tightened_to_equality
        parity_stats = parity_presolve!(model, propagator, postsolver, parity_stats_accumulator)
        model.infeasible && break
        parity_stats.domains_changed || break

        residue_stats = _residue_presolve_pass!(
            model,
            moduli,
            parity_stats.coefficient_changed_constraint_ids,
            postsolver,
            residue_stats_accumulator;
            treewidth_threshold = treewidth_threshold,
        )
    end

    return PresolveResult(
        model,
        postsolver,
        parity_stats_accumulator.stats,
        residue_stats_accumulator.stats,
    )
end
