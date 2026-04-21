const VarId = Int

"""
    StoredBitKind

Classify the bit fragments recorded during parity presolve.

Enumerate whether a stored bit is fixed, copied from a live binary variable, or
copied from its negation.
"""
@enum StoredBitKind FIXED0 FIXED1 BINVAR NEGATED_BINVAR

"""
    StoredBit

Store one reconstructed low-order bit.

Record the source of a bit appended during parity presolve, either as a fixed
value or as a reference to another binary variable.

# Fields
- `kind`: Bit source classification.
- `ref_var_id`: Referenced variable id for live binary bits, or `nothing` for
  fixed bits.
"""
struct StoredBit
    kind::StoredBitKind
    ref_var_id::Union{VarId, Nothing}
end

"""
    VarReconstruction

Store the reconstruction recipe for one original variable.

Track the surviving high-order variable, any fixed high-order fallback, an
additive offset, and the low-order bits accumulated during parity reductions.

# Fields
- `current_var_id`: Current reduced-model variable providing the high-order
  part, or `nothing` when that part is fixed.
- `fixed_high_order`: Fixed high-order value used when no live variable remains.
- `additive_offset`: Constant offset introduced by normalization or rewrites.
- `bits`: Stored low-order bits in least-significant-bit-first order.
"""
mutable struct VarReconstruction
    current_var_id::Union{VarId, Nothing}
    fixed_high_order::Union{Float64, Nothing}
    additive_offset::Float64
    bits::Vector{StoredBit}
end

"""
    ParityPostsolver

Store the data needed to reconstruct original variable values.

Keep reconstruction state for every tracked original variable so a solution of
the reduced parity-presolved model can be mapped back to the original space.

# Fields
- `original_var_ids`: Original variables to reconstruct, in sorted order.
- `var_data`: Reconstruction recipe for each tracked variable.

# See also
- [`postsolve`](@ref)
"""
struct ParityPostsolver
    original_var_ids::Vector{VarId}
    var_data::Dict{VarId, VarReconstruction}
end

"""
    _empty_reconstruction(var_id)

Create an empty reconstruction entry for `var_id`.

Initialize reconstruction state with `var_id` as the live high-order source,
zero additive offset, and no stored bits.
"""
@inline _empty_reconstruction(var_id::VarId) = VarReconstruction(var_id, nothing, 0.0, StoredBit[])

"""
    ParityPostsolver(var_ids)

Create a parity postsolver for a collection of original variables.

Sort and deduplicate `var_ids`, allocate reconstruction storage, and initialize
each tracked variable with an empty reconstruction recipe.

# Arguments
- `var_ids`: Iterable of variable ids that should be reconstructable after
  presolve.

# Returns
- A `ParityPostsolver` with one initialized entry per unique variable id.

# See also
- [`ensure_tracked_var!`](@ref)
- [`postsolve`](@ref)
"""
function ParityPostsolver(var_ids)
    original_var_ids = sort!(unique!(collect(var_ids)))
    var_data = Dict{VarId, VarReconstruction}()
    sizehint!(var_data, length(original_var_ids))

    for var_id in original_var_ids
        var_data[var_id] = _empty_reconstruction(var_id)
    end

    return ParityPostsolver(original_var_ids, var_data)
end

"""
    is_tracked_var(postsolver, var_id)

Check whether `var_id` is tracked by `postsolver`.

Return `true` when reconstruction data exists for `var_id`, and `false`
otherwise.
"""
is_tracked_var(postsolver::ParityPostsolver, var_id::VarId) = haskey(postsolver.var_data, var_id)

"""
    ensure_tracked_var!(postsolver, var_id)

Ensure that `postsolver` tracks `var_id`.

Insert an empty reconstruction entry when `var_id` is not already present.

# Arguments
- `postsolver`: Postsolver storage mutated in place.
- `var_id`: Variable id that must be tracked.

# Returns
- The mutated `postsolver`.

# Side Effects
- May add a new entry to `postsolver.var_data`.
"""
function ensure_tracked_var!(postsolver::ParityPostsolver, var_id::VarId)
    haskey(postsolver.var_data, var_id) || (postsolver.var_data[var_id] = _empty_reconstruction(var_id))
    return postsolver
end

"""
    _tracked_var_data!(postsolver, var_id)

Fetch mutable reconstruction data for `var_id`.

Ensure that `var_id` is tracked and return its reconstruction entry.

# Side Effects
- May add a new reconstruction entry to `postsolver`.
"""
function _tracked_var_data!(postsolver::ParityPostsolver, var_id::VarId)
    ensure_tracked_var!(postsolver, var_id)
    return postsolver.var_data[var_id]
end

"""
    _get_var_data(postsolver, var_id)

Return tracked reconstruction data for `var_id`.

Look up the reconstruction entry for `var_id` without creating one.

# Throws
- `ErrorException` if `var_id` is not tracked by `postsolver`.
"""
function _get_var_data(postsolver::ParityPostsolver, var_id::VarId)
    haskey(postsolver.var_data, var_id) || error("untracked postsolve variable $var_id")
    return postsolver.var_data[var_id]
end

"""
    append_fixed_bit!(postsolver, var_id, val)

Append a fixed low-order bit to the reconstruction of `var_id`.

Record `val` as the next stored bit in least-significant-bit-first order.

# Arguments
- `postsolver`: Postsolver storage mutated in place.
- `var_id`: Variable whose reconstruction receives the new bit.
- `val`: Fixed bit value to append.

# Returns
- The mutated `postsolver`.

# Side Effects
- Mutates the stored bit vector for `var_id`.
"""
function append_fixed_bit!(postsolver::ParityPostsolver, var_id::VarId, val::Bool)
    push!(_tracked_var_data!(postsolver, var_id).bits, StoredBit(val ? FIXED1 : FIXED0, nothing))
    return postsolver
end

"""
    append_binary_bit!(postsolver, var_id, ref_var_id; negated=false)

Append a live binary reference bit to the reconstruction of `var_id`.

Record that the next low-order bit of `var_id` should be read from
`ref_var_id`, optionally negated.

# Arguments
- `postsolver`: Postsolver storage mutated in place.
- `var_id`: Variable whose reconstruction receives the new bit.
- `ref_var_id`: Binary variable supplying the bit value at postsolve time.

# Keywords
- `negated=false`: Store `1 - ref_var_id` instead of `ref_var_id` when `true`.

# Returns
- The mutated `postsolver`.

# Side Effects
- Ensures both `var_id` and `ref_var_id` are tracked.
- Mutates the stored bit vector for `var_id`.
"""
function append_binary_bit!(
    postsolver::ParityPostsolver,
    var_id::VarId,
    ref_var_id::VarId;
    negated::Bool = false,
)
    ensure_tracked_var!(postsolver, ref_var_id)
    kind = negated ? NEGATED_BINVAR : BINVAR
    push!(_tracked_var_data!(postsolver, var_id).bits, StoredBit(kind, ref_var_id))
    return postsolver
end

"""
    mark_fixed_high_order!(postsolver, current_var_id, value)

Replace a live high-order variable with a fixed value.

Find every tracked reconstruction currently using `current_var_id` as its
high-order source and replace that source with the fixed value `value`.

# Arguments
- `postsolver`: Postsolver storage mutated in place.
- `current_var_id`: Live reduced-model variable being removed.
- `value`: Fixed high-order value to store instead.

# Returns
- The mutated `postsolver`.

# Side Effects
- Mutates matching reconstruction entries in `postsolver`.
"""
function mark_fixed_high_order!(postsolver::ParityPostsolver, current_var_id::VarId, value::Float64)
    for var_data in values(postsolver.var_data)
        var_data.current_var_id == current_var_id || continue
        var_data.current_var_id = nothing
        var_data.fixed_high_order = value
    end

    return postsolver
end

"""
    register_fixed_var!(postsolver, var_id, value)

Register that `var_id` has become fixed during presolve.

Alias [`mark_fixed_high_order!`](@ref) for call sites that speak in terms of a
fixed variable instead of a high-order source.
"""
register_fixed_var!(postsolver::ParityPostsolver, var_id::VarId, value::Float64) = mark_fixed_high_order!(postsolver, var_id, value)

"""
    add_reconstruction_offset!(postsolver, var_id, offset)

Add a constant reconstruction offset to `var_id`.

Accumulate `offset` into the constant term applied after rebuilding the
high-order and stored-bit components of `var_id`.

# Arguments
- `postsolver`: Postsolver storage mutated in place.
- `var_id`: Variable whose reconstruction offset should be updated.
- `offset`: Constant value to add.

# Returns
- The mutated `postsolver`.

# Side Effects
- Mutates the reconstruction entry for `var_id` when `offset != 0.0`.
"""
function add_reconstruction_offset!(postsolver::ParityPostsolver, var_id::VarId, offset::Float64)
    offset == 0.0 && return postsolver
    _tracked_var_data!(postsolver, var_id).additive_offset += offset
    return postsolver
end

"""
    _resolve_var_value(postsolver, var_id, solution, cache, active)

Resolve the reconstructed value of `var_id`.

Recursively evaluate the reconstruction recipe for `var_id`, using `solution`
for live reduced-model variables and `cache` to memoize completed values.

# Arguments
- `postsolver`: Reconstruction data source.
- `var_id`: Original or helper variable to reconstruct.
- `solution`: Reduced-model solution values keyed by variable id.
- `cache`: Memoization table for already reconstructed variables.
- `active`: Variables currently being resolved, used for cycle detection.

# Returns
- The reconstructed value of `var_id` as `Float64`.

# Throws
- `ErrorException` if a required reduced-model value is missing.
- `ErrorException` if reconstruction encounters a dependency cycle.

# Notes
- Variables absent from `postsolver.var_data` are treated as direct lookups in
  `solution`.
"""
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

"""
    _evaluate_stored_bit(bit, postsolver, solution, cache, active)

Evaluate one stored reconstruction bit.

Convert `bit` into its numeric contribution by resolving any referenced binary
variable and applying the stored bit-kind semantics.

# Returns
- The bit value as `0.0` or `1.0`.

# Throws
- `ErrorException` if a referenced variable id is missing.
- `ErrorException` if a referenced bit reconstructs to a non-binary value.
- `ErrorException` if `bit.kind` is unsupported.
"""
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

"""
    _high_order_value(postsolver, var_id, var_data, solution, cache, active)

Resolve the high-order component of a reconstructed variable.

Return the live reduced-model value referenced by `var_data.current_var_id`
when available, otherwise fall back to the fixed high-order value stored in
`var_data`.

# Returns
- The high-order contribution as `Float64`.

# Throws
- `ErrorException` if no live or fixed high-order value is available.
"""
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

"""
    postsolve(postsolver, solution)

Reconstruct the original-variable solution from a reduced parity solution.

Evaluate every tracked original variable in `postsolver` using the reduced-model
assignment `solution`, recursively resolving helper-variable dependencies and
stored low-order bits.

# Arguments
- `postsolver`: Reconstruction data produced during parity presolve.
- `solution`: Reduced-model solution values keyed by variable id.

# Returns
- A `Dict{VarId, Float64}` mapping each original variable id to its
  reconstructed value.

# Throws
- `ErrorException` if reconstruction requires a missing reduced-model value.
- `ErrorException` if reconstruction encounters inconsistent binary-bit data or
  dependency cycles.

# See also
- [`ParityPostsolver`](@ref)
"""
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
