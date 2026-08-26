Base.@kwdef mutable struct ParityStats
    parity_presolve_time::Float64 = 0.0
    num_parity_constraints_generated::Int = 0
    num_xor_constraints_generated::Int = 0
    num_xorand_constraints_generated::Int = 0
    num_parity_variables::Int = 0
    num_conjunction_variables::Int = 0
    num_total_boolean_variables::Int = 0
    num_parity_fixed_elimination::Int = 0
    num_parity_fixed_propagation::Int = 0
    num_parity_fixed_path::Int = 0
    num_constraints_split::Int = 0
    num_implication_edges_generated::Int = 0
    num_pattern_relations::Int = 0
    num_fixed_parity_substitutions::Int = 0
    num_pattern_substitutions::Int = 0
    num_parity_substitutions::Int = 0
    num_variables_fixed_after_parity_substitution::Int = 0
    num_constraints_modified_by_parity::Int = 0
    num_quadratic_terms_removed_by_parity::Int = 0
    num_constraints_linearized_by_parity::Int = 0
    num_constraints_removed_by_parity::Int = 0
    num_parity_presolve_rounds::Int = 0
    num_parity_presolve_phases::Int = 0
    infeasibility_detected::Bool = false
    infeasibility_source::Symbol = :none
end

mutable struct _ParityStatsAccumulator
    stats::ParityStats
    parity_variable_ids::Set{VarId}
    conjunction_variable_ids::Set{Tuple{VarId, VarId}}
    tracked_constraint_ids::Set{Int}
    modified_constraint_ids::Set{Int}
    linearized_constraint_ids::Set{Int}
    removed_constraint_ids::Set{Int}
    fixed_after_substitution_variable_ids::Set{VarId}
end

function _ParityStatsAccumulator(stats::ParityStats = ParityStats())
    return _ParityStatsAccumulator(
        stats,
        Set{VarId}(),
        Set{Tuple{VarId, VarId}}(),
        Set{Int}(),
        Set{Int}(),
        Set{Int}(),
        Set{Int}(),
        Set{VarId}(),
    )
end

_parity_stats(accumulator::_ParityStatsAccumulator) = accumulator.stats
_parity_stats(::Nothing) = ParityStats()

function _record_parity_presolve_time!(accumulator::_ParityStatsAccumulator, elapsed::Float64)
    accumulator.stats.parity_presolve_time += elapsed
    return accumulator
end

_record_parity_presolve_time!(::Nothing, ::Float64) = nothing

function _record_parity_presolve_round!(accumulator::_ParityStatsAccumulator)
    accumulator.stats.num_parity_presolve_rounds += 1
    return accumulator
end

_record_parity_presolve_round!(::Nothing) = nothing

function _record_parity_presolve_phase!(accumulator::_ParityStatsAccumulator)
    accumulator.stats.num_parity_presolve_phases += 1
    return accumulator
end

_record_parity_presolve_phase!(::Nothing) = nothing

function _sync_boolean_variable_counts!(accumulator::_ParityStatsAccumulator)
    stats = accumulator.stats
    stats.num_parity_variables = length(accumulator.parity_variable_ids)
    stats.num_conjunction_variables = length(accumulator.conjunction_variable_ids)
    stats.num_total_boolean_variables = stats.num_parity_variables + stats.num_conjunction_variables
    return accumulator
end

_sync_boolean_variable_counts!(::Nothing) = nothing

function _record_parity_variable!(accumulator::_ParityStatsAccumulator, var_id::VarId)
    push!(accumulator.parity_variable_ids, var_id)
    return _sync_boolean_variable_counts!(accumulator)
end

_record_parity_variable!(::Nothing, ::VarId) = nothing

function _record_conjunction_variable!(
    accumulator::_ParityStatsAccumulator,
    var_id_1::VarId,
    var_id_2::VarId,
)
    var_id_1 > var_id_2 && ((var_id_1, var_id_2) = (var_id_2, var_id_1))
    push!(accumulator.conjunction_variable_ids, (var_id_1, var_id_2))
    push!(accumulator.parity_variable_ids, var_id_1)
    push!(accumulator.parity_variable_ids, var_id_2)
    return _sync_boolean_variable_counts!(accumulator)
end

_record_conjunction_variable!(::Nothing, ::VarId, ::VarId) = nothing

function _record_generated_xor_constraint!(
    accumulator::_ParityStatsAccumulator,
    support::AbstractVector{Bool},
    pos_to_var_id::AbstractVector{VarId},
)
    stats = accumulator.stats
    stats.num_parity_constraints_generated += 1
    stats.num_xor_constraints_generated += 1
    for idx in eachindex(support)
        support[idx] || continue
        push!(accumulator.parity_variable_ids, pos_to_var_id[idx])
    end
    return _sync_boolean_variable_counts!(accumulator)
end

_record_generated_xor_constraint!(::Nothing, ::AbstractVector{Bool}, ::AbstractVector{VarId}) = nothing

function _record_builder_generated!(accumulator::_ParityStatsAccumulator, builder::XorConstraintBuilder)
    stats = accumulator.stats
    stats.num_parity_constraints_generated += 1

    has_conjunction = false
    for (var_id, is_active) in builder.par
        is_active || continue
        push!(accumulator.parity_variable_ids, var_id)
    end
    for ((var_id_1, var_id_2), is_active) in builder.conj
        is_active || continue
        has_conjunction = true
        _record_conjunction_variable!(accumulator, var_id_1, var_id_2)
    end

    if has_conjunction
        stats.num_xorand_constraints_generated += 1
    else
        stats.num_xor_constraints_generated += 1
    end

    return _sync_boolean_variable_counts!(accumulator)
end

_record_builder_generated!(::Nothing, ::XorConstraintBuilder) = nothing

function _record_implication_edge!(accumulator::_ParityStatsAccumulator, added::Bool)
    added && (accumulator.stats.num_implication_edges_generated += 1)
    return accumulator
end

_record_implication_edge!(::Nothing, ::Bool) = nothing

function _record_parity_fixing!(accumulator::_ParityStatsAccumulator, source::Symbol)
    if source == :elimination
        accumulator.stats.num_parity_fixed_elimination += 1
    elseif source == :propagation
        accumulator.stats.num_parity_fixed_propagation += 1
    elseif source == :path
        accumulator.stats.num_parity_fixed_path += 1
    end
    return accumulator
end

_record_parity_fixing!(::Nothing, ::Symbol) = nothing

function _record_constraint_split!(accumulator::_ParityStatsAccumulator)
    accumulator.stats.num_constraints_split += 1
    return accumulator
end

_record_constraint_split!(::Nothing) = nothing

function _record_fixed_parity_substitution!(accumulator::_ParityStatsAccumulator)
    stats = accumulator.stats
    stats.num_fixed_parity_substitutions += 1
    stats.num_parity_substitutions += 1
    return accumulator
end

_record_fixed_parity_substitution!(::Nothing) = nothing

function _record_pattern_substitution!(
    accumulator::_ParityStatsAccumulator,
    component_size::Int,
    rewritten_variables::Int,
)
    stats = accumulator.stats
    stats.num_pattern_relations += max(0, component_size - 1)
    stats.num_pattern_substitutions += rewritten_variables
    stats.num_parity_substitutions += rewritten_variables
    return accumulator
end

_record_pattern_substitution!(::Nothing, ::Int, ::Int) = nothing

function _record_infeasibility!(accumulator::_ParityStatsAccumulator, source::Symbol)
    stats = accumulator.stats
    stats.infeasibility_detected = true
    stats.infeasibility_source == :none && (stats.infeasibility_source = source)
    return accumulator
end

_record_infeasibility!(::Nothing, ::Symbol) = nothing
