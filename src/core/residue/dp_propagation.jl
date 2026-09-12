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

struct _StandardizedResidueConstraint
    con::Constraint
    var_bounds::Dict{VarId, IntVar}
    var_ubs::Dict{VarId, Int}
    constraint_shift::Int
    objective_shift::Float64
end

struct _ResidueSetResult
    residues::BitVector
    saturated::Bool
end

struct _ResidueCacheEntry
    modulus::Int
    residues::BitVector
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
        ;
        treewidth_threshold::Integer = typemax(Int),
    )
    td = minimum_degree_tree_decomposition(component)
    if tree_decomposition_width(td) > treewidth_threshold
        return ResidueDPResult(component, trues(modulus), true)
    end

    stack = ConditionedResidueSet[]
    removed_ids = Set{VarId}()

    for action in action_order(component, td)
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
        ;
        treewidth_threshold::Integer = typemax(Int),
    )
    m > 0 || throw(ArgumentError("modulus must be positive, got $m"))
    modulus = Int(m)
    capped_var_ubs = _validate_residue_bounds(modulus, con, var_bounds)
    interaction_graph = InteractionGraph(con, modulus)

    results = ResidueDPResult[]
    for component in decompose(interaction_graph)
        component isa NonSingleton || continue
        push!(
            results,
            _component_residue_dp(
                modulus,
                component,
                capped_var_ubs;
                treewidth_threshold = treewidth_threshold,
            ),
        )
    end
    return results
end

function _validate_residue_domain(var_id::VarId, var::IntVar)
    isfinite(var.lb) ||
        throw(ArgumentError("variable $var_id must have finite lower bound, got $(var.lb)"))
    isfinite(var.ub) ||
        throw(ArgumentError("variable $var_id must have finite upper bound, got $(var.ub)"))
    isinteger(var.lb) ||
        throw(ArgumentError("variable $var_id must have integer lower bound, got $(var.lb)"))
    isinteger(var.ub) ||
        throw(ArgumentError("variable $var_id must have integer upper bound, got $(var.ub)"))
    var.lb <= var.ub ||
        throw(ArgumentError("variable $var_id has inconsistent bounds [$(var.lb), $(var.ub)]"))
    return nothing
end

function _has_supported_residue_domain(var::IntVar)
    return isfinite(var.lb) &&
        isfinite(var.ub) &&
        isinteger(var.lb) &&
        isinteger(var.ub)
end

function _standardize_residue_constraint(
        con::Constraint,
        var_bounds::Dict{VarId, IntVar},
        obj_expr::Union{Nothing, QuadExpr} = nothing,
    )
    is_integer(con) ||
        throw(ArgumentError("residue propagation requires an integer-valued constraint"))

    shifted_con = deepcopy(con)
    shift_expr = deepcopy(con.qe)
    shifted_obj = obj_expr === nothing ? nothing : deepcopy(obj_expr)
    shifted_bounds = Dict{VarId, IntVar}()
    shifted_ubs = Dict{VarId, Int}()
    var_ids = sort!(collect(vars(con.qe)))

    for var_id in var_ids
        haskey(var_bounds, var_id) ||
            throw(ArgumentError("missing bounds for variable $var_id"))
        var = var_bounds[var_id]
        _validate_residue_domain(var_id, var)

        lb = trunc(Int, var.lb)
        ub = trunc(Int, var.ub)
        width = ub - lb
        offset = Float64(lb)

        affine_transform!(shifted_con, var_id, 1.0, offset)
        affine_transform!(shift_expr, var_id, 1.0, offset)
        shifted_obj !== nothing && affine_transform!(shifted_obj, var_id, 1.0, offset)

        shifted_bounds[var_id] = IntVar(0.0, Float64(width))
        shifted_ubs[var_id] = width
    end

    isinteger(shift_expr.constant) ||
        throw(ArgumentError("standardized constraint shift must be integer-valued, got $(shift_expr.constant)"))
    constraint_shift = trunc(Int, shift_expr.constant)
    objective_shift = shifted_obj === nothing ? 0.0 : shifted_obj.constant - obj_expr.constant

    return _StandardizedResidueConstraint(
        shifted_con,
        shifted_bounds,
        shifted_ubs,
        constraint_shift,
        objective_shift,
    )
end

function _generate_residue_moduli(strategy::Symbol, threshold::Integer)
    threshold_int = Int(threshold)

    if strategy == :small_primes
        threshold_int < 2 && return Int[]
        return [candidate for candidate in 2:threshold_int if _is_prime(candidate)]
    elseif strategy == :powers_of_two
        threshold_int < 2 && return Int[]
        moduli = Int[]
        modulus = 2
        while modulus <= threshold_int
            push!(moduli, modulus)
            modulus > div(typemax(Int), 2) && break
            modulus *= 2
        end
        return moduli
    elseif strategy == :divisor_free
        threshold_int < 2 && return Int[]
        moduli = Int[]
        for candidate in threshold_int:-1:2
            _is_divisor_free_candidate(candidate, moduli) || continue
            push!(moduli, candidate)
        end
        return sort!(moduli)
    end

    throw(ArgumentError("unsupported residue modulus strategy $strategy"))
end

function _is_divisor_free_candidate(candidate::Int, moduli::AbstractVector{Int})
    for modulus in moduli
        (modulus % candidate == 0 || candidate % modulus == 0) && return false
    end
    return true
end

function _is_prime(candidate::Int)
    candidate < 2 && return false
    candidate == 2 && return true
    iseven(candidate) && return false

    divisor = 3
    while divisor <= div(candidate, divisor)
        candidate % divisor == 0 && return false
        divisor += 2
    end
    return true
end

function _zero_residue_set(modulus::Int)
    residues = falses(modulus)
    residues[1] = true
    return residues
end

function _convolve_linear_singleton!(
        residues::BitVector,
        component::LinSingleton,
        ub::Integer,
        modulus::Int,
    )
    capped_ub = min(Int(ub), modulus - 1)
    capped_ub >= 0 ||
        throw(ArgumentError("linear singleton upper bound must be nonnegative, got $ub"))

    nvalues = capped_ub + 1
    span = 1
    step = mod(component.lin_coeff, modulus)

    while 2 * span <= nvalues
        _or_shifted_residues!(residues, copy(residues), span * step, modulus)
        _is_saturated(residues) && return residues
        span *= 2
    end

    remaining = nvalues - span
    if remaining > 0
        _or_shifted_residues!(residues, copy(residues), remaining * step, modulus)
    end

    return residues
end

function _quad_singleton_residues(
        component::QuadSingleton,
        ub::Integer,
        modulus::Int,
    )
    capped_ub = min(Int(ub), modulus - 1)
    capped_ub >= 0 ||
        throw(ArgumentError("quadratic singleton upper bound must be nonnegative, got $ub"))

    residues = falses(modulus)
    @inbounds for value in 0:capped_ub
        residue = mod(
            component.quad_coeff * value * value + component.lin_coeff * value,
            modulus,
        )
        residues[residue + 1] = true
    end
    return residues
end

function _convolve_residue_set(
        residues::BitVector,
        component_residues::BitVector,
        modulus::Int,
    )
    _is_saturated(component_residues) && return trues(modulus)
    return _modular_sumset(residues, component_residues, modulus)
end

function _compute_achievable_residues(
        m::Integer,
        con::Constraint,
        var_bounds::Dict{VarId, IntVar},
        ;
        treewidth_threshold::Integer = typemax(Int),
    )
    m > 0 || throw(ArgumentError("modulus must be positive, got $m"))
    modulus = Int(m)
    capped_var_ubs = _validate_residue_bounds(modulus, con, var_bounds)
    components = decompose(InteractionGraph(con, modulus))
    residues = _zero_residue_set(modulus)

    for component in components
        component isa LinSingleton || continue
        _convolve_linear_singleton!(
            residues,
            component,
            capped_var_ubs[component.var_id],
            modulus,
        )
        _is_saturated(residues) && return _ResidueSetResult(trues(modulus), true)
    end

    for component in components
        component isa QuadSingleton || continue
        component_residues = _quad_singleton_residues(
            component,
            capped_var_ubs[component.var_id],
            modulus,
        )
        residues = _convolve_residue_set(residues, component_residues, modulus)
        _is_saturated(residues) && return _ResidueSetResult(trues(modulus), true)
    end

    for component in components
        component isa NonSingleton || continue
        result = _component_residue_dp(
            modulus,
            component,
            capped_var_ubs;
            treewidth_threshold = treewidth_threshold,
        )
        residues = _convolve_residue_set(residues, result.residues, modulus)
        _is_saturated(residues) && return _ResidueSetResult(trues(modulus), true)
    end

    return _ResidueSetResult(residues, _is_saturated(residues))
end

function _smallest_attainable_ge(
        lower::Int,
        residues::BitVector,
        modulus::Int,
        constraint_shift::Int,
    )
    @inbounds for delta in 0:(modulus - 1)
        candidate = lower + delta
        residues[mod(candidate - constraint_shift, modulus) + 1] && return candidate
    end
    error("nonempty residue set has no attainable lower bound modulo $modulus")
end

function _largest_attainable_le(
        upper::Int,
        residues::BitVector,
        modulus::Int,
        constraint_shift::Int,
    )
    @inbounds for delta in 0:(modulus - 1)
        candidate = upper - delta
        residues[mod(candidate - constraint_shift, modulus) + 1] && return candidate
    end
    error("nonempty residue set has no attainable upper bound modulo $modulus")
end

function _tighten_constraint_bounds_by_residues!(
        con::Constraint,
        residues::BitVector,
        modulus::Int,
        constraint_shift::Int,
    )
    changed = false

    if isfinite(con.lhs)
        new_lhs = Float64(_smallest_attainable_ge(
            ceil(Int, con.lhs),
            residues,
            modulus,
            constraint_shift,
        ))
        if new_lhs > con.lhs
            con.lhs = _canonicalize_zero(new_lhs)
            changed = true
        end
    end

    if isfinite(con.rhs)
        new_rhs = Float64(_largest_attainable_le(
            floor(Int, con.rhs),
            residues,
            modulus,
            constraint_shift,
        ))
        if new_rhs < con.rhs
            con.rhs = _canonicalize_zero(new_rhs)
            changed = true
        end
    end

    return changed
end

function _reapply_residue_cache!(
        con::Constraint,
        cache::Vector{_ResidueCacheEntry},
        constraint_shift::Int,
    )
    changed = false
    while true
        pass_changed = false
        for entry in cache
            pass_changed = _tighten_constraint_bounds_by_residues!(
                con,
                entry.residues,
                entry.modulus,
                constraint_shift,
            ) || pass_changed
            con.lhs > con.rhs && return true
        end
        changed = changed || pass_changed
        pass_changed || return changed
    end
end

function _residue_constraint_status(
        con::Constraint,
        var_bounds::Dict{VarId, IntVar},
    )
    con.lhs > con.rhs && return :infeasible
    is_integer(con) || return :skip

    for var_id in vars(con.qe)
        haskey(var_bounds, var_id) ||
            throw(ArgumentError("missing bounds for variable $var_id"))
        var = var_bounds[var_id]
        var.lb > var.ub && return :infeasible
        _has_supported_residue_domain(var) || return :skip
    end

    return :process
end

function _residue_constraint_result(
        changed::Bool,
        tightened_to_equality::Bool,
        processed::Bool,
        infeasible::Bool,
        stats_accumulator::Union{Nothing, _ResidueStatsAccumulator},
    )
    return (
        changed = changed,
        tightened_to_equality = tightened_to_equality,
        processed = processed,
        infeasible = infeasible,
        residue_stats = _residue_stats(stats_accumulator),
    )
end

"""
    residue_presolve_constraint!(model, con, moduli; treewidth_threshold=typemax(Int), collect_stats=false)
    residue_presolve_constraint!(model, con, moduli, stats_accumulator; treewidth_threshold=typemax(Int))

Run residue-based bound tightening on a single constraint of `model`.

# Returns
- A named tuple with `changed`, `tightened_to_equality`, `processed`, and
  `infeasible` fields, plus accumulated `residue_stats`.

# Side Effects
- May tighten `con.lhs` or `con.rhs`.
- May set `model.infeasible = true`.
"""
function residue_presolve_constraint!(
        model::QPModel,
        con::Constraint,
        moduli::AbstractVector{<:Integer};
        treewidth_threshold::Integer = typemax(Int),
        collect_stats::Bool = false,
    )
    stats_accumulator = collect_stats ? _ResidueStatsAccumulator() : nothing
    return residue_presolve_constraint!(
        model,
        con,
        moduli,
        stats_accumulator;
        treewidth_threshold = treewidth_threshold,
    )
end

function residue_presolve_constraint!(
        model::QPModel,
        con::Constraint,
        moduli::AbstractVector{<:Integer},
        stats_accumulator::Union{Nothing, _ResidueStatsAccumulator};
        treewidth_threshold::Integer = typemax(Int),
    )
    start_time = time()
    try
        return _residue_presolve_constraint_impl!(
            model,
            con,
            moduli,
            stats_accumulator;
            treewidth_threshold = treewidth_threshold,
        )
    finally
        _record_residue_presolve_time!(stats_accumulator, time() - start_time)
    end
end

function _residue_presolve_constraint_impl!(
        model::QPModel,
        con::Constraint,
        moduli::AbstractVector{<:Integer},
        stats_accumulator::Union{Nothing, _ResidueStatsAccumulator};
        treewidth_threshold::Integer,
    )
    model.infeasible && return _residue_constraint_result(
        false,
        false,
        false,
        true,
        stats_accumulator,
    )

    was_inequality = !is_equality(con)
    status = _residue_constraint_status(con, model.vars)
    if status == :infeasible
        model.infeasible = true
        return _residue_constraint_result(
            false,
            false,
            false,
            true,
            stats_accumulator,
        )
    elseif status == :skip
        return _residue_constraint_result(
            false,
            false,
            false,
            false,
            stats_accumulator,
        )
    end
    isempty(moduli) && return _residue_constraint_result(
        false,
        false,
        false,
        false,
        stats_accumulator,
    )

    standardized = _standardize_residue_constraint(con, model.vars, model.obj_expr)
    cache = _ResidueCacheEntry[]
    changed = false
    processed = false

    for modulus_value in moduli
        modulus = Int(modulus_value)
        result = _compute_achievable_residues(
            modulus,
            standardized.con,
            standardized.var_bounds;
            treewidth_threshold = treewidth_threshold,
        )
        _record_residue_modulus_evaluated!(stats_accumulator, con)
        processed = true
        result.saturated && continue

        push!(cache, _ResidueCacheEntry(modulus, result.residues))
        old_lhs = con.lhs
        old_rhs = con.rhs
        was_equality_before_reapply = is_equality(con)
        reapply_changed = _reapply_residue_cache!(con, cache, standardized.constraint_shift)
        bounds_changed = con.lhs != old_lhs || con.rhs != old_rhs
        bounds_changed && _record_residue_bound_tightening!(
            stats_accumulator,
            con,
            old_lhs,
            old_rhs,
            was_equality_before_reapply,
        )
        residue_infeasible = con.lhs > con.rhs
        (bounds_changed || residue_infeasible) &&
            _record_useful_residue_modulus!(stats_accumulator)
        changed = reapply_changed || bounds_changed || changed

        if residue_infeasible
            model.infeasible = true
            _record_residue_infeasibility!(stats_accumulator)
            return _residue_constraint_result(
                changed,
                false,
                processed,
                true,
                stats_accumulator,
            )
        end
    end

    return _residue_constraint_result(
        changed,
        changed && was_inequality && is_equality(con),
        processed,
        false,
        stats_accumulator,
    )
end

"""
    residue_presolve!(model, strategy; threshold, treewidth_threshold=typemax(Int)) -> QPModel

Run residue-based constraint-bound propagation on `model`.

The pass first normalizes `model`, then applies residue reductions to the
remaining constraints. Nonlinear components whose tree-decomposition width
exceeds `treewidth_threshold` are treated as saturated for that modulus.
"""
function residue_presolve!(
        model::QPModel,
        strategy::Symbol;
        threshold::Integer,
        treewidth_threshold::Integer = typemax(Int),
    )::QPModel
    moduli = _generate_residue_moduli(strategy, threshold)
    _residue_presolve_pass!(
        model,
        moduli,
        nothing;
        treewidth_threshold = treewidth_threshold,
    )
    return model
end
