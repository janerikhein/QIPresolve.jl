"""
    ConditionedResidueSet

Store achievable residues modulo `modulus`, conditioned on assignments to
`var_ids`.

The assignment index is lexicographic over sorted `var_ids`: for assignment
`x`, `index = 1 + sum(x[i] * strides[i])`, with the last variable changing
fastest.
"""
mutable struct ConditionedResidueSet
    modulus::Int
    var_ids::Vector{VarId}
    var_ubs::Vector{Int}
    strides::Vector{Int}
    residues::Vector{BitVector}
    local_contribution::Vector{Int}
end

"""
    ResidueDPResult

Result of the nonlinear-component residue dynamic program.
"""
struct ResidueDPResult
    component::NonSingleton
    residues::BitVector
    saturated::Bool
end

function _assignment_count(var_ubs::AbstractVector{Int})
    count = 1
    for ub in var_ubs
        count *= ub + 1
    end
    return count
end

function _compute_strides(var_ubs::AbstractVector{Int})
    strides = Vector{Int}(undef, length(var_ubs))
    next_stride = 1
    for index in length(var_ubs):-1:1
        strides[index] = next_stride
        next_stride *= var_ubs[index] + 1
    end
    return strides
end

function _assignment_to_index(
        assignment::AbstractVector{Int},
        strides::AbstractVector{Int},
    )
    index = 1
    @inbounds for i in eachindex(assignment, strides)
        index += assignment[i] * strides[i]
    end
    return index
end

function _index_to_assignment!(
        assignment::Vector{Int},
        index::Int,
        var_ubs::AbstractVector{Int},
        strides::AbstractVector{Int},
    )
    offset = index - 1
    @inbounds for i in eachindex(var_ubs, strides)
        assignment[i] = fld(offset, strides[i])
        offset -= assignment[i] * strides[i]
    end
    return assignment
end

function _copy_without(values::Vector{T}, skipped::Int) where {T}
    output = Vector{T}(undef, length(values) - 1)
    output_index = 1
    @inbounds for index in eachindex(values)
        index == skipped && continue
        output[output_index] = values[index]
        output_index += 1
    end
    return output
end

function _is_saturated(bits::BitVector)
    return all(bits)
end

function _has_saturated_residue_set(set::ConditionedResidueSet)
    return any(_is_saturated, set.residues)
end

function _or_shifted_residues!(
        destination::BitVector,
        source::BitVector,
        shift::Int,
        modulus::Int,
    )
    shift_mod = mod(shift, modulus)
    @inbounds for residue in 0:(modulus - 1)
        source[residue + 1] || continue
        destination[mod(residue + shift_mod, modulus) + 1] = true
    end
    return destination
end

function _modular_sumset(left::BitVector, right::BitVector, modulus::Int)
    output = falses(modulus)
    @inbounds for left_residue in 0:(modulus - 1)
        left[left_residue + 1] || continue
        for right_residue in 0:(modulus - 1)
            right[right_residue + 1] || continue
            output[mod(left_residue + right_residue, modulus) + 1] = true
        end
    end
    return output
end

function ConditionedResidueSet(modulus::Integer)
    modulus > 0 || throw(ArgumentError("modulus must be positive, got $modulus"))
    m = Int(modulus)
    residues = [falses(m)]
    residues[1][1] = true
    return ConditionedResidueSet(
        m,
        VarId[],
        Int[],
        Int[],
        residues,
        zeros(Int, 1),
    )
end

"""
    condition_on!(set, id, ub)

Add `id` to the conditioning variables with domain `{0, ..., ub}`.
"""
function condition_on!(set::ConditionedResidueSet, id::VarId, ub::Integer)
    ub_int = Int(ub)
    0 <= ub_int <= set.modulus - 1 ||
        throw(ArgumentError("upper bound for $id must be in 0:$(set.modulus - 1), got $ub"))

    insert_position = searchsortedfirst(set.var_ids, id)
    if insert_position <= length(set.var_ids) && set.var_ids[insert_position] == id
        throw(ArgumentError("variable $id is already conditioned on"))
    end

    old_var_ids = set.var_ids
    old_var_ubs = set.var_ubs
    old_strides = set.strides
    old_residues = set.residues

    new_var_ids = Vector{VarId}(undef, length(old_var_ids) + 1)
    new_var_ubs = Vector{Int}(undef, length(old_var_ubs) + 1)
    @inbounds for new_position in eachindex(new_var_ids)
        if new_position < insert_position
            new_var_ids[new_position] = old_var_ids[new_position]
            new_var_ubs[new_position] = old_var_ubs[new_position]
        elseif new_position == insert_position
            new_var_ids[new_position] = id
            new_var_ubs[new_position] = ub_int
        else
            new_var_ids[new_position] = old_var_ids[new_position - 1]
            new_var_ubs[new_position] = old_var_ubs[new_position - 1]
        end
    end

    new_strides = _compute_strides(new_var_ubs)
    new_residues = Vector{BitVector}(undef, _assignment_count(new_var_ubs))

    old_assignment = zeros(Int, length(old_var_ids))
    for old_index in eachindex(old_residues)
        _index_to_assignment!(old_assignment, old_index, old_var_ubs, old_strides)

        base_index = 1
        old_position = 1
        @inbounds for new_position in eachindex(new_var_ids)
            if new_position == insert_position
                continue
            end
            base_index += old_assignment[old_position] * new_strides[new_position]
            old_position += 1
        end

        @inbounds for value in 0:ub_int
            new_index = base_index + value * new_strides[insert_position]
            new_residues[new_index] = copy(old_residues[old_index])
        end
    end

    set.var_ids = new_var_ids
    set.var_ubs = new_var_ubs
    set.strides = new_strides
    set.residues = new_residues
    set.local_contribution = zeros(Int, length(new_residues))
    return set
end

"""
    populate_local_contribution!(set, target_id, quad_terms, lin_terms, removed_ids)

Populate local contributions for the remove action of `target_id`.
"""
function populate_local_contribution!(
        set::ConditionedResidueSet,
        target_id::VarId,
        quad_terms,
        lin_terms,
        removed_ids,
    )
    target_position = findfirst(==(target_id), set.var_ids)
    target_position === nothing &&
        throw(ArgumentError("target variable $target_id is not conditioned on"))

    modulus = set.modulus
    var_positions = Dict(var_id => position for (position, var_id) in enumerate(set.var_ids))
    diagonal = mod(get(quad_terms, (target_id, target_id), 0), modulus)
    linear = mod(get(lin_terms, target_id, 0), modulus)
    neighbor_coeffs = zeros(Int, length(set.var_ids))

    for ((first_id, second_id), coefficient) in quad_terms
        first_id == second_id && continue
        if first_id == target_id
            other_id = second_id
        elseif second_id == target_id
            other_id = first_id
        else
            continue
        end

        other_id in removed_ids && continue
        other_position = get(var_positions, other_id, 0)
        other_position == 0 &&
            throw(ArgumentError("active neighbor $other_id of $target_id is not conditioned on"))
        neighbor_coeffs[other_position] =
            mod(neighbor_coeffs[other_position] + coefficient, modulus)
    end

    local_contribution = Vector{Int}(undef, length(set.residues))
    values = zeros(Int, length(set.var_ids))
    directions = ones(Int, length(set.var_ids))
    index = 1
    target_value = 0
    neighbor_sum = 0
    contribution = 0

    @inbounds for iteration in eachindex(local_contribution)
        local_contribution[index] = contribution
        iteration == length(local_contribution) && break

        changed_position = 0
        delta = 0
        for position in length(values):-1:1
            next_value = values[position] + directions[position]
            if 0 <= next_value <= set.var_ubs[position]
                changed_position = position
                delta = directions[position]
                values[position] = next_value
                index += delta * set.strides[position]
                break
            end
            directions[position] = -directions[position]
        end
        changed_position == 0 && error("failed to advance mixed-radix Gray code")

        if changed_position == target_position
            old_target_value = target_value
            target_value += delta
            contribution = mod(
                contribution +
                diagonal * (target_value * target_value - old_target_value * old_target_value) +
                linear * delta +
                delta * neighbor_sum,
                modulus,
            )
        else
            coefficient = neighbor_coeffs[changed_position]
            if coefficient != 0
                neighbor_sum = mod(neighbor_sum + coefficient * delta, modulus)
                contribution = mod(
                    contribution + target_value * coefficient * delta,
                    modulus,
                )
            end
        end
    end

    set.local_contribution = local_contribution
    return set
end

"""
    fold_on!(set, target_id)

Shift residues by local contribution and remove `target_id` by OR-merging over
its domain.
"""
function fold_on!(set::ConditionedResidueSet, target_id::VarId)
    target_position = findfirst(==(target_id), set.var_ids)
    target_position === nothing &&
        throw(ArgumentError("target variable $target_id is not conditioned on"))
    length(set.local_contribution) == length(set.residues) ||
        throw(ArgumentError("local contributions must be populated before folding"))

    new_var_ids = _copy_without(set.var_ids, target_position)
    new_var_ubs = _copy_without(set.var_ubs, target_position)
    new_strides = _compute_strides(new_var_ubs)
    new_residues = [falses(set.modulus) for _ in 1:_assignment_count(new_var_ubs)]

    old_assignment = zeros(Int, length(set.var_ids))
    for old_index in eachindex(set.residues)
        _index_to_assignment!(old_assignment, old_index, set.var_ubs, set.strides)

        new_index = 1
        new_position = 1
        @inbounds for old_position in eachindex(set.var_ids)
            old_position == target_position && continue
            new_index += old_assignment[old_position] * new_strides[new_position]
            new_position += 1
        end

        _or_shifted_residues!(
            new_residues[new_index],
            set.residues[old_index],
            set.local_contribution[old_index],
            set.modulus,
        )
    end

    set.var_ids = new_var_ids
    set.var_ubs = new_var_ubs
    set.strides = new_strides
    set.residues = new_residues
    set.local_contribution = zeros(Int, length(new_residues))
    return set
end

"""
    join!(left, right)

Join two conditioned residue sets by modular sumset for each conditioning
assignment.
"""
function join!(left::ConditionedResidueSet, right::ConditionedResidueSet)
    left.modulus == right.modulus ||
        throw(ArgumentError("cannot join residue sets with different moduli"))
    left.var_ids == right.var_ids ||
        throw(ArgumentError("cannot join residue sets with different conditioned variables"))
    left.var_ubs == right.var_ubs ||
        throw(ArgumentError("cannot join residue sets with different variable domains"))
    length(left.residues) == length(right.residues) ||
        throw(ArgumentError("cannot join residue sets with different assignment counts"))

    modulus = left.modulus
    for index in eachindex(left.residues, right.residues)
        left.residues[index] = _modular_sumset(
            left.residues[index],
            right.residues[index],
            modulus,
        )
    end
    left.local_contribution = zeros(Int, length(left.residues))
    return left
end

condition_on(set::ConditionedResidueSet, id::VarId, ub::Integer) =
    condition_on!(set, id, ub)
populate_local_contribution(
    set::ConditionedResidueSet,
    target_id::VarId,
    quad_terms,
    lin_terms,
    removed_ids,
) = populate_local_contribution!(
    set,
    target_id,
    quad_terms,
    lin_terms,
    removed_ids,
)
fold_on(set::ConditionedResidueSet, target_id::VarId) = fold_on!(set, target_id)

function _validate_residue_bounds(
        modulus::Int,
        con::Constraint,
        var_bounds::Dict{VarId, IntVar},
    )
    var_ubs = Dict{VarId, Int}()
    for var_id in vars(con.qe)
        haskey(var_bounds, var_id) ||
            throw(ArgumentError("missing bounds for variable $var_id"))
        var = var_bounds[var_id]
        var.lb == 0.0 ||
            throw(ArgumentError("variable $var_id must have lower bound 0, got $(var.lb)"))
        isfinite(var.ub) ||
            throw(ArgumentError("variable $var_id must have finite upper bound, got $(var.ub)"))
        isinteger(var.ub) ||
            throw(ArgumentError("variable $var_id must have integer upper bound, got $(var.ub)"))
        var.ub >= 0.0 ||
            throw(ArgumentError("variable $var_id must have nonnegative upper bound, got $(var.ub)"))
        var_ubs[var_id] = min(trunc(Int, var.ub), modulus - 1)
    end
    return var_ubs
end

function _component_residue_dp(
        modulus::Int,
        component::NonSingleton,
        var_ubs::Dict{VarId, Int},
    )
    stack = ConditionedResidueSet[]
    removed_ids = Set{VarId}()

    for action in action_order(component)
        if action isa LeafAction
            push!(stack, ConditionedResidueSet(modulus))
        elseif action isa IntroduceAction
            isempty(stack) && error("introduce action requires an existing residue set")
            condition_on!(last(stack), action.id, var_ubs[action.id])
        elseif action isa RemoveAction
            isempty(stack) && error("remove action requires an existing residue set")
            state = last(stack)
            populate_local_contribution!(
                state,
                action.id,
                component.quad_coeffs,
                component.lin_coeffs,
                removed_ids,
            )
            fold_on!(state, action.id)
            push!(removed_ids, action.id)
        elseif action isa JoinAction
            length(stack) >= 2 || error("join action requires two residue sets")
            right = pop!(stack)
            left = pop!(stack)
            join!(left, right)
            push!(stack, left)
        else
            error("unknown residue action $(typeof(action))")
        end

        if _has_saturated_residue_set(last(stack))
            return ResidueDPResult(component, trues(modulus), true)
        end
    end

    length(stack) == 1 || error("residue DP ended with $(length(stack)) states")
    final_state = only(stack)
    isempty(final_state.var_ids) || error("residue DP root state must be unconditioned")
    final_residues = copy(only(final_state.residues))
    return ResidueDPResult(component, final_residues, _is_saturated(final_residues))
end

"""
    compute_nonlinear_residue_sets(m, con, var_bounds) -> Vector{ResidueDPResult}

Compute achievable residue sets modulo `m` for each nonlinear interaction
component of `con`.
"""
function compute_nonlinear_residue_sets(
        m::Integer,
        con::Constraint,
        var_bounds::Dict{VarId, IntVar},
    )
    m > 0 || throw(ArgumentError("modulus must be positive, got $m"))
    modulus = Int(m)
    capped_var_ubs = _validate_residue_bounds(modulus, con, var_bounds)
    interaction_graph = InteractionGraph(con, modulus)

    results = ResidueDPResult[]
    for component in decompose(interaction_graph)
        component isa NonSingleton || continue
        push!(results, _component_residue_dp(modulus, component, capped_var_ubs))
    end
    return results
end
