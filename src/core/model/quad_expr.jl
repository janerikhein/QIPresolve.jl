const VarId = Int

"""
    QuadExpr

Semi-sparse representation of a quadratic expression

    x'Qx + c'x + constant

with support for incremental variable addition variable removal.

# Design

`QuadExpr` stores quadratic and linear coefficients in *dense buffers* of
size `cap × cap` and `cap`, but only the first `nvars` logical variables are
active at any time.

To avoid expensive row/column copying on removals, the structure separates
**logical variable positions** from **physical buffer indices**:

- Logical positions `1:nvars` represent the active variables in the expression
- Physical buffer indices `1:cap` index into `quad_buf` and `lin_buf`
- `perm[pos]` maps a logical position to its physical buffer index

This indirection allows variable removal by swapping logical positions and
shrinking `nvars` without modifying the underlying buffers.

# Fields

- `nvars::Int`  
  Number of active variables in the expression.

- `cap::Int`  
  Capacity of the internal buffers. Physical slots `nvars+1:cap` are reserved
  for future variable additions.

- `pos_to_var::Vector{VarId}`  
  Maps logical position `pos ∈ 1:nvars` to the corresponding variable id.

- `var_to_pos::Dict{VarId,Int}`  
  Inverse mapping from variable id to logical position.

- `perm::Vector{Int}`  
  Lazy permutation mapping logical positions to physical buffer indices.

- `quad_buf::Matrix{Float64}`  
  Dense `cap × cap` buffer storing quadratic coefficients, indexed by physical
  buffer indices.

- `lin_buf::Vector{Float64}`  
  Dense buffer of length `cap` storing linear coefficients, indexed by physical
  buffer indices.

- `constant::Float64`  
  Constant term of the quadratic expression.

# Access patterns

Quadratic and linear coefficients are accessed via:

    pos_i = var_to_pos[id_i]
    phys_i = perm[pos_i]

The active quadratic and linear parts can be viewed as dense objects via:

    quad(qe) :: AbstractMatrix{Float64}
    lin(qe)  :: AbstractVector{Float64}

which return views into the underlying buffers.

# Operations

- Variable addition and removal are O(1) (amortized for addition due to geometric buffer growth)
  with optional buffer clearing in O(n). Default behaviour is to clear buffer on removal but not on addition.
- Coefficient access and modification are O(1).

This structure is intended for presolve and model-transformation workflows
where expressions are modified frequently and structural changes must be cheap.
"""
mutable struct QuadExpr
    nvars::Int
    cap::Int
    pos_to_var::Vector{VarId}
    var_to_pos::Dict{VarId, Int}
    perm::Vector{Int}
    quad_buf::Matrix{Float64}
    lin_buf::Vector{Float64}
    constant::Float64
end

"""
    QuadExpr(quad_terms, lin_terms; constant=0.0, buf_size_factor=2.0) -> QuadExpr

Construct a `QuadExpr` from lists of quadratic and linear terms.

`quad_terms` is a vector of `(coef, id1, id2)` tuples, and `lin_terms` is a
vector of `(coef, id)` tuples. Duplicate variables are deduplicated, and
duplicate terms are accumulated.

Buffer capacity is set to `ceil(buf_size_factor * nvars)` (at least 1) to allow
future variable additions without immediate resizing.
"""
function QuadExpr(
        quad_terms::Vector{Tuple{Float64, VarId, VarId}},
        lin_terms::Vector{Tuple{Float64, VarId}};
        constant::Float64 = 0.0,
        buf_size_factor::Float64 = 2.0,
    )
    @assert buf_size_factor >= 1.0

    # collect variables
    vars = Vector{VarId}()
    seen = Set{VarId}()

    for (_, i, j) in quad_terms
        if !(i in seen)
            push!(vars, i)
            push!(seen, i)
        end
        if !(j in seen)
            push!(vars, j)
            push!(seen, j)
        end
    end
    for (_, i) in lin_terms
        if !(i in seen)
            push!(vars, i); push!(seen, i)
        end
    end

    nvars = length(vars)
    cap = max(1, Int(ceil(buf_size_factor * nvars)))

    # Allocate buffers and mappings
    pos_to_var = Vector{VarId}(undef, cap)
    var_to_pos = Dict{VarId, Int}()
    perm = collect(1:cap)

    quad_buf = zeros(Float64, cap, cap)
    lin_buf = zeros(Float64, cap)

    # register variables at logical positions 1:nvars
    for (pos, v) in enumerate(vars)
        pos_to_var[pos] = v
        var_to_pos[v] = pos
    end

    qe = QuadExpr(
        nvars,
        cap,
        pos_to_var,
        var_to_pos,
        perm,
        quad_buf,
        lin_buf,
        constant,
    )

    # insert quadratic terms
    for (coef, id1, id2) in quad_terms
        pos1 = qe.var_to_pos[id1]
        pos2 = qe.var_to_pos[id2]
        i = qe.perm[pos1]
        j = qe.perm[pos2]
        if i <= j
            @inbounds qe.quad_buf[i, j] += coef
        else
            @inbounds qe.quad_buf[j, i] += coef
        end
    end

    # insert linear terms
    for (coef, id) in lin_terms
        pos = qe.var_to_pos[id]
        i = qe.perm[pos]
        @inbounds qe.lin_buf[i] += coef
    end

    return qe
end

"""
    quad(qe) -> AbstractMatrix{Float64}

Return a view of the active quadratic coefficient matrix for `qe`.

The returned view indexes the first `nvars` logical variables via the current
`perm` mapping, so it stays consistent under variable additions/removals.
Mutating the view writes into the underlying `quad_buf`.
"""
@inline function quad(qe::QuadExpr)
    p = @view qe.perm[1:qe.nvars]
    return @view qe.quad_buf[p, p]
end

"""
    lin(qe) -> AbstractVector{Float64}

Return a view of the active linear coefficient vector for `qe`.

The returned view indexes the first `nvars` logical variables via the current
`perm` mapping, so it stays consistent under variable additions/removals.
Mutating the view writes into the underlying `lin_buf`.
"""
@inline function lin(qe::QuadExpr)
    p = @view qe.perm[1:qe.nvars]
    return @view qe.lin_buf[p]
end


"""
    checks if expression is empty, assumes expr is normalized
"""
is_empty(qe::QuadExpr) = (qe.nvars == 0)

"""
    checks if expression is single linear term, i.e. var bound, assumes expr is normalized
"""
is_singleton(qe::QuadExpr) = (qe.nvars == 1 && qe.quad_buf[qe.perm[1], qe.perm[1]] == 0)

"""
    vars(qe) -> AbstractVector{VarId}

Return a view of all variables currently present in the expression.
The returned object is allocation-free and reflects the current state.
"""
@inline function vars(qe::QuadExpr)
    return @view qe.pos_to_var[1:qe.nvars]
end

"""
    eval_full(qe, x) -> Float64

Evaluate the full quadratic expression for the dense vector `x` indexed by `VarId`.
"""
function eval_full(qe::QuadExpr, x::AbstractVector{<:Real})
    total = qe.constant
    for posi in 1:qe.nvars
        idi = qe.pos_to_var[posi]
        xi = x[idi]
        pi = qe.perm[posi]
        total += qe.lin_buf[pi] * xi
        for posj in 1:qe.nvars
            idj = qe.pos_to_var[posj]
            xj = x[idj]
            pj = qe.perm[posj]
            total += qe.quad_buf[pi, pj] * xi * xj
        end
    end
    return total
end

"""
    has_var(qe, id) -> Bool

Return `true` if variable `id` is present in the quadratic expression,
`false` otherwise.
"""
@inline function has_var(qe::QuadExpr, id::VarId)::Bool
    return haskey(qe.var_to_pos, id)
end


"""
    var_pos(qe, id) -> Int

Return the logical position of variable `id` in `qe`.

Returns `0` if the variable is not present.
"""
@inline function var_pos(qe::QuadExpr, id::VarId)::Int
    return get(qe.var_to_pos, id, 0)  # 0 means "not present"
end


"""
    get_quad_coeff(qe, id1, id2) -> Float64

Return the quadratic coefficient for variables `id1` and `id2`.

If either variable is not present in the expression, returns `0.0`.
"""
@inline function get_quad_coeff(qe::QuadExpr, id1::VarId, id2::VarId)
    pos1 = var_pos(qe, id1)
    pos1 == 0 && return 0.0
    pos2 = var_pos(qe, id2)
    pos2 == 0 && return 0.0

    i = qe.perm[pos1]
    j = qe.perm[pos2]
    if i <= j
        @inbounds return qe.quad_buf[i, j]
    else
        @inbounds return qe.quad_buf[j, i]
    end
end

"""
    get_lin_coeff(qe, id)

Return the linear coefficient for variable `id`.

If the variable is not present in the expression, returns `0.0`.
"""
@inline function get_lin_coeff(qe::QuadExpr, id::VarId)
    pos = var_pos(qe, id)
    pos == 0 && return 0.0
    i = qe.perm[pos]
    @inbounds return qe.lin_buf[i]
end


"""
    set_lin_coeff!(qe, id, val) -> Bool

Set the linear coefficient for variable `id` to `val`.

Returns `true` if `id` exists in the expression, otherwise returns `false`
(and does not modify the expression).

Note: Updates `lin_buf` at the permuted physical index.
"""
function set_lin_coeff!(qe::QuadExpr, id::VarId, val::Float64)::Bool
    pos = var_pos(qe, id)
    pos == 0 && return false
    i = qe.perm[pos]  # physical
    @inbounds qe.lin_buf[i] = val
    return true
end

"""
    add_lin_coeff!(qe, id, delta) -> Bool

Add `delta` to the linear coefficient of variable `id`.

Returns `true` if `id` exists in the expression, otherwise returns `false`
(and does not modify the expression).
"""
function add_lin_coeff!(qe::QuadExpr, id::VarId, delta::Float64)::Bool
    pos = var_pos(qe, id)
    pos == 0 && return false
    i = qe.perm[pos]  # physical
    @inbounds qe.lin_buf[i] += delta
    return true
end


"""
    set_quad_coeff!(qe, id1, id2, val) -> Bool

Set the quadratic coefficient for variables `id1` and `id2` to `val`.

Returns `true` if both variables exist in the expression, otherwise returns `false`
(and does not modify the expression).

Note: Coefficients are stored canonically in the upper triangle of `quad_buf`
(by physical buffer index). Lower-triangular entries are forced to zero.
The `sym` keyword is accepted for compatibility but has no effect.
"""
function set_quad_coeff!(qe::QuadExpr, id1::VarId, id2::VarId, val::Float64; sym::Bool = false)
    pos1 = var_pos(qe, id1)
    pos1 == 0 && return false
    pos2 = var_pos(qe, id2)
    pos2 == 0 && return false

    i = qe.perm[pos1]  # physical
    j = qe.perm[pos2]  # physical

    if i <= j
        @inbounds qe.quad_buf[i, j] = val
        i != j && (@inbounds qe.quad_buf[j, i] = 0.0)
    else
        @inbounds qe.quad_buf[j, i] = val
        @inbounds qe.quad_buf[i, j] = 0.0
    end
    return true
end

"""
    add_quad_coeff!(qe, id1, id2, val) -> Bool

Set the quadratic coefficient for variables `id1` and `id2` to `val`.

Returns `true` if both variables exist in the expression, otherwise returns `false`
(and does not modify the expression).

Note: Coefficients are stored canonically in the upper triangle of `quad_buf`
(by physical buffer index). Lower-triangular entries are forced to zero.
The `sym` keyword is accepted for compatibility but has no effect.
"""
function add_quad_coeff!(qe::QuadExpr, id1::VarId, id2::VarId, delta::Float64; sym::Bool = false)
    pos1 = var_pos(qe, id1)
    pos1 == 0 && return false
    pos2 = var_pos(qe, id2)
    pos2 == 0 && return false

    i = qe.perm[pos1]  # physical
    j = qe.perm[pos2]  # physical

    if i <= j
        @inbounds qe.quad_buf[i, j] += delta
        i != j && (@inbounds qe.quad_buf[j, i] = 0.0)
    else
        @inbounds qe.quad_buf[j, i] += delta
        @inbounds qe.quad_buf[i, j] = 0.0
    end
    return true
end


"""
    remove_var!(qe, id; clear_buf=true) -> Bool

Remove variable `id` from the expression.

Returns `true` if `id` exists, otherwise returns `false` and leaves `qe` unchanged.

When `clear_buf=true` (default), zeroes the freed physical row/column and linear
slot in the buffers. This is O(n) in the buffer size but keeps old coefficients
from being accidentally reused if the physical slot is later reactivated.
"""
function remove_var!(qe::QuadExpr, id::VarId; clear_buf::Bool = true)
    pos = var_pos(qe, id)
    pos == 0 && return false

    last = qe.nvars
    phys_pos = qe.perm[pos]
    phys_last = qe.perm[last]

    if pos != last
        last_id = qe.pos_to_var[last]
        qe.pos_to_var[pos] = last_id
        qe.var_to_pos[last_id] = pos

        qe.perm[pos] = phys_last
    end

    # put freed physical slot into inactive tail (so add_var! can reuse it)
    qe.perm[last] = phys_pos

    delete!(qe.var_to_pos, id)

    # optional clearing buffer
    if clear_buf
        @inbounds begin
            qe.lin_buf[phys_pos] = 0.0
            @views begin
                fill!(qe.quad_buf[phys_pos, :], 0.0)
                fill!(qe.quad_buf[:, phys_pos], 0.0)
            end
        end
    end

    qe.nvars -= 1
    return true
end

function normalize!(qe::QuadExpr)
    for pos in reverse(1:qe.nvars)
        phys_pos = qe.perm[pos]

        qe.lin_buf[phys_pos] != 0 && continue

        all_zero = true
        for other_pos in 1:qe.nvars
            phys_other = qe.perm[other_pos]
            q = if phys_pos <= phys_other
                @inbounds qe.quad_buf[phys_pos, phys_other]
            else
                @inbounds qe.quad_buf[phys_other, phys_pos]
            end
            if q != 0
                all_zero = false
                break
            end
        end
        if all_zero
            remove_var!(qe, qe.pos_to_var[pos])
        end
    end
    return
end


"""
    _ensure_capacity!(qe, needed)

Ensure internal buffers can hold at least `needed` logical variables.

Grows `quad_buf`, `lin_buf`, and the position mappings geometrically, preserving
existing content and initializing new inactive logical positions to map to their
own physical slots.
"""
function _ensure_capacity!(qe::QuadExpr, needed::Int)
    needed <= qe.cap && return

    oldcap = qe.cap
    newcap = oldcap
    while newcap < needed
        newcap *= 2
    end

    # grow buffers
    new_quad = zeros(Float64, newcap, newcap)
    new_lin = zeros(Float64, newcap)
    new_pos_to_var = Vector{VarId}(undef, newcap)
    new_perm = Vector{Int}(undef, newcap)

    # copy old content
    @inbounds begin
        new_quad[1:oldcap, 1:oldcap] .= qe.quad_buf
        new_lin[1:oldcap] .= qe.lin_buf
        new_pos_to_var[1:oldcap] .= qe.pos_to_var
        new_perm[1:oldcap] .= qe.perm
    end

    # initialize new (inactive) logical positions to map to new physical slots
    @inbounds for k in (oldcap + 1):newcap
        new_perm[k] = k
        # new_pos_to_var[k] left undef because it's inactive
    end

    qe.quad_buf = new_quad
    qe.lin_buf = new_lin
    qe.pos_to_var = new_pos_to_var
    qe.perm = new_perm
    qe.cap = newcap
    return
end

"""
    add_var!(qe, id) -> Int

Add variable `id` if it is not present.
Returns the logical position (1:qe.nvars) of the variable.

Uses a free physical slot from `perm[qe.nvars+1]` and increments `nvars`.
Grows capacity geometrically when needed.
"""
function add_var!(qe::QuadExpr, id::VarId; clear_buf::Bool = false)
    pos = var_pos(qe, id)
    pos != 0 && return pos

    _ensure_capacity!(qe, qe.nvars + 1)

    newpos = qe.nvars + 1
    phys = qe.perm[newpos]   # free physical slot

    # register variable at logical position newpos
    qe.nvars = newpos
    qe.pos_to_var[newpos] = id
    qe.var_to_pos[id] = newpos

    # optional clearing buffer
    if clear_buf
        @inbounds begin
            qe.lin_buf[phys] = 0.0
            @views begin
                fill!(qe.quad_buf[phys, :], 0.0)
                fill!(qe.quad_buf[:, phys], 0.0)
            end
        end
    end

    return newpos
end

"""
    invert_affine(scale::Float64, offset::Float64) -> (inv_scale, inv_offset)

Return the affine parameters that undo `x -> scale * x + offset`.

Throws `ArgumentError` when `scale == 0.0`.
"""
@inline function invert_affine(scale::Float64, offset::Float64)
    scale == 0.0 && throw(ArgumentError("Cannot invert affine transform with scale = 0"))
    inv_scale = 1.0 / scale
    inv_offset = -offset / scale
    return inv_scale, inv_offset
end

"""
    affine_transform!(
        qe::QuadExpr,
        var_id::VarId,
        scale::Float64,
        offset::Float64;
        invert::Bool = false,
    ) -> QuadExpr

Apply the affine substitution to one variable:

    x_var_id := scale * x_var_id + offset

for the quadratic expression represented by `qe`:

    ∑_{i,j} Q[i,j] x_i x_j + ∑_i c[i] x_i + constant

Notes
- This updates `quad_buf`, `lin_buf`, and `constant` in-place.
- Assumes canonical upper-triangular quadratic storage.
- No allocations; does not create any views.
- If `invert=true`, applies the inverse transform for the provided `(scale, offset)`.
"""
function affine_transform!(qe::QuadExpr, var_id::VarId, scale::Float64, offset::Float64; invert::Bool = false)
    if invert
        scale, offset = invert_affine(scale, offset)
    end

    posk = var_pos(qe, var_id)
    posk == 0 && return qe

    pk = qe.perm[posk]  # physical index of var_id

    @inbounds begin
        ck = qe.lin_buf[pk]
        qkk = qe.quad_buf[pk, pk]

        # constant update from ck*xk + qkk*xk^2 with xk := scale*xk + offset
        qe.constant += offset * ck + (offset * offset) * qkk

        # update linear coefficient of xk
        qe.lin_buf[pk] = scale * ck + 2.0 * scale * offset * qkk

        # update interactions with all other variables
        for posj in 1:qe.nvars
            pj = qe.perm[posj]
            pj == pk && continue

            qkj = if pk <= pj
                @inbounds qe.quad_buf[pk, pj]
            else
                @inbounds qe.quad_buf[pj, pk]
            end

            if qkj != 0.0
                qe.lin_buf[pj] += offset * qkj
                if pk <= pj
                    @inbounds qe.quad_buf[pk, pj] = scale * qkj
                    pk != pj && (@inbounds qe.quad_buf[pj, pk] = 0.0)
                else
                    @inbounds qe.quad_buf[pj, pk] = scale * qkj
                    @inbounds qe.quad_buf[pk, pj] = 0.0
                end
            end
        end

        # update diagonal entry
        qe.quad_buf[pk, pk] = (scale * scale) * qkk
    end

    return qe
end

"""
    invert_lin(a::Float64, b::Float64) -> (a_inv, b_inv)

Return the linear parameters that undo `x_k -> a * x_k + b * x_m`.

Throws `ArgumentError` when `a == 0.0`.
"""
@inline function invert_lin(a::Float64, b::Float64)
    a == 0.0 && throw(ArgumentError("Cannot invert shear x <- a*x + b*y with a = 0"))
    a_inv = 1.0 / a
    b_inv = -b / a
    return a_inv, b_inv
end

@inline function invert_lin(a::Float64, b::Float64, c::Float64)
    a == 0.0 && throw(ArgumentError("Cannot invert affine shear x <- a*x + b*y + c with a = 0"))
    a_inv = 1.0 / a
    b_inv = -b / a
    c_inv = -c / a
    return a_inv, b_inv, c_inv
end


"""
    lin_transform!(
        qe::QuadExpr,
        var_id::VarId,
        other_id::VarId,
        a::Float64,
        b::Float64;
        c::Float64 = 0.0,
        invert::Bool = false,
    ) -> QuadExpr

Perform the linear 2-variable substitution

    x_var_id := a * x_var_id + b * x_other_id + c

in the quadratic expression

    ∑_{i,j} Q[i,j] x_i x_j + ∑_i c[i] x_i + constant

stored in `qe`.

Notes
- Updates `quad_buf` and `lin_buf` in-place (no allocations, no views).
- `other_id` (the "y") is not otherwise changed; only `var_id` is substituted.
- Returns unchanged if `var_id` is absent.
- Treats self-substitution (`var_id == other_id`) as a no-op.
- Requires `other_id` to already exist unless `b == 0.0`, in which case the
  substitution reduces to an affine transform on `var_id`.
- If `invert=true`, applies the inverse transform for the provided `(a, b, c)`.
- Coefficients are stored canonically in the upper triangle of `quad_buf`.
"""
function lin_transform!(
        qe::QuadExpr,
        var_id::VarId,
        other_id::VarId,
        a::Float64,
        b::Float64;
        c::Float64 = 0.0,
        invert::Bool = false
    )
    if invert
        (a, b, c) = invert_lin(a, b, c)
    end

    posk = var_pos(qe, var_id)
    posk == 0 && return qe
    if var_id == other_id
        return qe
    end

    posm = var_pos(qe, other_id)
    if posm == 0
        return b == 0.0 ? affine_transform!(qe, var_id, a, c) : qe
    end

    pk = qe.perm[posk]  # physical index of var_id (k)
    pm = qe.perm[posm]  # physical index of other_id (m)

    pk == pm && return qe

    @inbounds begin
        # Save original coefficients involving k and m before overwriting.
        ck = qe.lin_buf[pk]
        cm = qe.lin_buf[pm]

        qkk = qe.quad_buf[pk, pk]
        qkm = if pk <= pm
            @inbounds qe.quad_buf[pk, pm]
        else
            @inbounds qe.quad_buf[pm, pk]
        end
        qmm = qe.quad_buf[pm, pm]

        qe.constant += c * ck + (c * c) * qkk

        # Linear terms
        qe.lin_buf[pk] = a * ck + 2.0 * a * c * qkk
        qe.lin_buf[pm] = cm + b * ck + c * qkm + 2.0 * b * c * qkk

        # Quadratic terms among {k,m}
        qe.quad_buf[pk, pk] = (a * a) * qkk
        qkmpm = a * qkm + (2.0 * a * b) * qkk
        if pk <= pm
            @inbounds qe.quad_buf[pk, pm] = qkmpm
            pk != pm && (@inbounds qe.quad_buf[pm, pk] = 0.0)
        else
            @inbounds qe.quad_buf[pm, pk] = qkmpm
            @inbounds qe.quad_buf[pk, pm] = 0.0
        end
        qe.quad_buf[pm, pm] = qmm + b * qkm + (b * b) * qkk

        # Interactions with all other vars j ≠ k,m
        for posj in 1:qe.nvars
            pj = qe.perm[posj]
            (pj == pk || pj == pm) && continue

            qkj_old = if pk <= pj
                @inbounds qe.quad_buf[pk, pj]
            else
                @inbounds qe.quad_buf[pj, pk]
            end

            if qkj_old != 0.0
                qmj_old = if pm <= pj
                    @inbounds qe.quad_buf[pm, pj]
                else
                    @inbounds qe.quad_buf[pj, pm]
                end

                qe.lin_buf[pj] += c * qkj_old

                qmj_new = qmj_old + b * qkj_old
                if pm <= pj
                    @inbounds qe.quad_buf[pm, pj] = qmj_new
                    pm != pj && (@inbounds qe.quad_buf[pj, pm] = 0.0)
                else
                    @inbounds qe.quad_buf[pj, pm] = qmj_new
                    @inbounds qe.quad_buf[pm, pj] = 0.0
                end

                qkj_new = a * qkj_old
                if pk <= pj
                    @inbounds qe.quad_buf[pk, pj] = qkj_new
                    pk != pj && (@inbounds qe.quad_buf[pj, pk] = 0.0)
                else
                    @inbounds qe.quad_buf[pj, pk] = qkj_new
                    @inbounds qe.quad_buf[pk, pj] = 0.0
                end
            end
        end
    end

    return qe
end

function scale_transform!(qe::QuadExpr, scale::Float64)
    for posi in 1:qe.nvars
        pi = qe.perm[posi]
        for posj in posi:qe.nvars
            pj = qe.perm[posj]
            q = if pi <= pj
                @inbounds qe.quad_buf[pi, pj]
            else
                @inbounds qe.quad_buf[pj, pi]
            end
            qscaled = scale * q
            if pi <= pj
                @inbounds qe.quad_buf[pi, pj] = qscaled
                pi != pj && (@inbounds qe.quad_buf[pj, pi] = 0.0)
            else
                @inbounds qe.quad_buf[pj, pi] = qscaled
                @inbounds qe.quad_buf[pi, pj] = 0.0
            end
        end
    end
    for pos in 1:qe.nvars
        i = qe.perm[pos]
        qe.lin_buf[i] *= scale
    end

    return qe
end
