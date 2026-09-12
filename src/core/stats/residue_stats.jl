Base.@kwdef mutable struct ResidueStats
    residue_presolve_time::Float64 = 0.0
    num_constraints_processed::Int = 0
    num_constraints_tightened::Int = 0
    num_constraints_eq_tightened::Int = 0
    num_lower_bounds_tightened::Int = 0
    num_upper_bounds_tightened::Int = 0
    infeasibility_detected::Bool = false
    avg_relative_bound_tightening::Float64 = 0.0
    num_moduli_evaluated::Int = 0
    num_useful_moduli::Int = 0
end

mutable struct _ResidueStatsAccumulator
    stats::ResidueStats
    processed_constraint_ids::Set{Int}
    tightened_constraint_ids::Set{Int}
    eq_tightened_constraint_ids::Set{Int}
    lower_bound_tightened_constraint_ids::Set{Int}
    upper_bound_tightened_constraint_ids::Set{Int}
    baseline_bounds::Dict{Int, Tuple{Float64, Float64}}
    relative_tightening_by_constraint_id::Dict{Int, Float64}
end

function _ResidueStatsAccumulator(stats::ResidueStats = ResidueStats())
    return _ResidueStatsAccumulator(
        stats,
        Set{Int}(),
        Set{Int}(),
        Set{Int}(),
        Set{Int}(),
        Set{Int}(),
        Dict{Int, Tuple{Float64, Float64}}(),
        Dict{Int, Float64}(),
    )
end

_residue_stats(accumulator::_ResidueStatsAccumulator) = accumulator.stats
_residue_stats(::Nothing) = ResidueStats()

function _record_residue_presolve_time!(accumulator::_ResidueStatsAccumulator, elapsed::Float64)
    accumulator.stats.residue_presolve_time += elapsed
    return accumulator
end

_record_residue_presolve_time!(::Nothing, ::Float64) = nothing

function _insert_new!(ids::Set{Int}, id::Int)
    old_length = length(ids)
    push!(ids, id)
    return length(ids) != old_length
end

function _ensure_residue_baseline!(accumulator::_ResidueStatsAccumulator, con::Constraint)
    haskey(accumulator.baseline_bounds, con.id) ||
        (accumulator.baseline_bounds[con.id] = (con.lhs, con.rhs))
    return accumulator
end

_ensure_residue_baseline!(::Nothing, ::Constraint) = nothing

function _sync_relative_bound_tightening!(accumulator::_ResidueStatsAccumulator)
    values_iter = values(accumulator.relative_tightening_by_constraint_id)
    accumulator.stats.avg_relative_bound_tightening = isempty(values_iter) ?
        0.0 :
        sum(values_iter; init = 0.0) / length(values_iter)
    return accumulator
end

_sync_relative_bound_tightening!(::Nothing) = nothing

function _record_residue_processed!(accumulator::_ResidueStatsAccumulator, con::Constraint)
    _ensure_residue_baseline!(accumulator, con)
    _insert_new!(accumulator.processed_constraint_ids, con.id) &&
        (accumulator.stats.num_constraints_processed += 1)
    return accumulator
end

_record_residue_processed!(::Nothing, ::Constraint) = nothing

function _record_residue_modulus_evaluated!(accumulator::_ResidueStatsAccumulator, con::Constraint)
    _record_residue_processed!(accumulator, con)
    accumulator.stats.num_moduli_evaluated += 1
    return accumulator
end

_record_residue_modulus_evaluated!(::Nothing, ::Constraint) = nothing

function _record_useful_residue_modulus!(accumulator::_ResidueStatsAccumulator)
    accumulator.stats.num_useful_moduli += 1
    return accumulator
end

_record_useful_residue_modulus!(::Nothing) = nothing

function _record_residue_infeasibility!(accumulator::_ResidueStatsAccumulator)
    accumulator.stats.infeasibility_detected = true
    return accumulator
end

_record_residue_infeasibility!(::Nothing) = nothing

function _update_relative_bound_tightening!(
        accumulator::_ResidueStatsAccumulator,
        con::Constraint,
    )
    haskey(accumulator.baseline_bounds, con.id) || return accumulator
    baseline_lhs, baseline_rhs = accumulator.baseline_bounds[con.id]
    isfinite(baseline_lhs) && isfinite(baseline_rhs) || return accumulator

    baseline_range = baseline_rhs - baseline_lhs
    baseline_range > 0.0 || return accumulator

    improvement = 0.0
    isfinite(con.lhs) && (improvement += max(0.0, con.lhs - baseline_lhs))
    isfinite(con.rhs) && (improvement += max(0.0, baseline_rhs - con.rhs))
    improvement > 0.0 || return accumulator

    accumulator.relative_tightening_by_constraint_id[con.id] = improvement / baseline_range
    return _sync_relative_bound_tightening!(accumulator)
end

_update_relative_bound_tightening!(::Nothing, ::Constraint) = nothing

function _record_residue_bound_tightening!(
        accumulator::_ResidueStatsAccumulator,
        con::Constraint,
        old_lhs::Float64,
        old_rhs::Float64,
        was_equality::Bool,
    )
    lower_tightened = isfinite(old_lhs) && con.lhs > old_lhs
    upper_tightened = isfinite(old_rhs) && con.rhs < old_rhs
    (lower_tightened || upper_tightened) || return accumulator

    _insert_new!(accumulator.tightened_constraint_ids, con.id) &&
        (accumulator.stats.num_constraints_tightened += 1)
    lower_tightened &&
        _insert_new!(accumulator.lower_bound_tightened_constraint_ids, con.id) &&
        (accumulator.stats.num_lower_bounds_tightened += 1)
    upper_tightened &&
        _insert_new!(accumulator.upper_bound_tightened_constraint_ids, con.id) &&
        (accumulator.stats.num_upper_bounds_tightened += 1)

    if !was_equality && is_equality(con)
        _insert_new!(accumulator.eq_tightened_constraint_ids, con.id) &&
            (accumulator.stats.num_constraints_eq_tightened += 1)
    end

    return _update_relative_bound_tightening!(accumulator, con)
end

_record_residue_bound_tightening!(::Nothing, ::Constraint, ::Float64, ::Float64, ::Bool) =
    nothing
