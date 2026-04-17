const VarId = Int

@enum StoredBitKind FIXED0 FIXED1 BINVAR NEGATED_BINVAR

struct StoredBit
    kind::StoredBitKind
    ref_var_id::Union{VarId, Nothing}
end

mutable struct VarReconstruction
    current_var_id::Union{VarId, Nothing}
    fixed_high_order::Union{Float64, Nothing}
    additive_offset::Float64
    bits::Vector{StoredBit}
end

struct ParityPostsolver
    original_var_ids::Vector{VarId}
    var_data::Dict{VarId, VarReconstruction}
end

function ParityPostsolver(var_ids)
    original_var_ids = sort!(unique!(collect(var_ids)))
    var_data = Dict{VarId, VarReconstruction}()
    sizehint!(var_data, length(original_var_ids))

    for var_id in original_var_ids
        var_data[var_id] = VarReconstruction(var_id, nothing, 0.0, StoredBit[])
    end

    return ParityPostsolver(original_var_ids, var_data)
end

is_tracked_var(postsolver::ParityPostsolver, var_id::VarId) = haskey(postsolver.var_data, var_id)

function ensure_tracked_var!(postsolver::ParityPostsolver, var_id::VarId)
    haskey(postsolver.var_data, var_id) || (postsolver.var_data[var_id] = VarReconstruction(var_id, nothing, 0.0, StoredBit[]))
    return postsolver
end

function _get_var_data(postsolver::ParityPostsolver, var_id::VarId)
    haskey(postsolver.var_data, var_id) || error("untracked postsolve variable $var_id")
    return postsolver.var_data[var_id]
end

function append_fixed_bit!(postsolver::ParityPostsolver, var_id::VarId, val::Bool)
    ensure_tracked_var!(postsolver, var_id)
    bit = StoredBit(val ? FIXED1 : FIXED0, nothing)
    push!(_get_var_data(postsolver, var_id).bits, bit)
    return postsolver
end

function append_binary_bit!(
    postsolver::ParityPostsolver,
    var_id::VarId,
    ref_var_id::VarId;
    negated::Bool = false,
)
    ensure_tracked_var!(postsolver, var_id)
    ensure_tracked_var!(postsolver, ref_var_id)
    kind = negated ? NEGATED_BINVAR : BINVAR
    push!(_get_var_data(postsolver, var_id).bits, StoredBit(kind, ref_var_id))
    return postsolver
end

function mark_fixed_high_order!(postsolver::ParityPostsolver, current_var_id::VarId, value::Float64)
    for var_data in values(postsolver.var_data)
        var_data.current_var_id == current_var_id || continue
        var_data.current_var_id = nothing
        var_data.fixed_high_order = value
    end

    return postsolver
end

function register_fixed_var!(postsolver::ParityPostsolver, var_id::VarId, value::Float64)
    mark_fixed_high_order!(postsolver, var_id, value)
    return postsolver
end

function add_reconstruction_offset!(postsolver::ParityPostsolver, var_id::VarId, offset::Float64)
    offset == 0.0 && return postsolver
    ensure_tracked_var!(postsolver, var_id)
    _get_var_data(postsolver, var_id).additive_offset += offset
    return postsolver
end

function _resolve_var_value(
    postsolver::ParityPostsolver,
    var_id::VarId,
    solution::AbstractDict{VarId, <:Real},
    cache::Dict{VarId, Float64},
    active::Set{VarId},
)
    haskey(cache, var_id) && return cache[var_id]
    haskey(postsolver.var_data, var_id) || begin
        haskey(solution, var_id) || error("missing reduced-solution value for variable $var_id")
        return Float64(solution[var_id])
    end
    var_id in active && error("cycle detected while reconstructing variable $var_id")

    push!(active, var_id)
    var_data = _get_var_data(postsolver, var_id)
    high_order = _high_order_value(postsolver, var_id, var_data, solution, cache, active)

    bit_offset = 0.0
    bit_weight = 1.0
    for bit in var_data.bits
        bit_offset += bit_weight * _evaluate_stored_bit(bit, postsolver, solution, cache, active)
        bit_weight *= 2.0
    end

    value = high_order * bit_weight + bit_offset + var_data.additive_offset
    cache[var_id] = value
    delete!(active, var_id)
    return value
end

function _evaluate_stored_bit(
    bit::StoredBit,
    postsolver::ParityPostsolver,
    solution::AbstractDict{VarId, <:Real},
    cache::Dict{VarId, Float64},
    active::Set{VarId},
)
    if bit.kind == FIXED0
        return 0.0
    elseif bit.kind == FIXED1
        return 1.0
    end

    ref_var_id = bit.ref_var_id
    ref_var_id === nothing && error("missing referenced variable id for stored bit")
    value = _resolve_var_value(postsolver, ref_var_id, solution, cache, active)
    (value == 0.0 || value == 1.0) || error("stored bit reference $ref_var_id reconstructed to non-binary value $value")

    if bit.kind == BINVAR
        return value
    elseif bit.kind == NEGATED_BINVAR
        return 1.0 - value
    end

    error("unsupported stored bit kind $(bit.kind)")
end

function _high_order_value(
    postsolver::ParityPostsolver,
    var_id::VarId,
    var_data::VarReconstruction,
    solution::AbstractDict{VarId, <:Real},
    cache::Dict{VarId, Float64},
    active::Set{VarId},
)
    current_var_id = var_data.current_var_id
    if current_var_id !== nothing
        if haskey(solution, current_var_id)
            return Float64(solution[current_var_id])
        elseif current_var_id != var_id && haskey(postsolver.var_data, current_var_id)
            return _resolve_var_value(postsolver, current_var_id, solution, cache, active)
        end
    end

    fixed_high_order = var_data.fixed_high_order
    fixed_high_order !== nothing || error("missing high-order value for reconstructed variable")
    return fixed_high_order
end

function postsolve(postsolver::ParityPostsolver, solution::AbstractDict{VarId, <:Real})
    original_solution = Dict{VarId, Float64}()
    sizehint!(original_solution, length(postsolver.original_var_ids))
    cache = Dict{VarId, Float64}()
    active = Set{VarId}()

    for original_var_id in postsolver.original_var_ids
        original_solution[original_var_id] = _resolve_var_value(postsolver, original_var_id, solution, cache, active)
    end

    return original_solution
end
