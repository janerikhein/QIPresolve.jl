"""
    XorConstraintMeta

Store cached metadata for an `XorConstraint`.

# Fields
- `nnz_par`: Number of active parity terms.
- `nnz_conj`: Number of active conjunctive terms, counted once per unordered
  pair.
- `is_pure_xor`: Whether the constraint has no conjunctive terms.
- `requires_update`: Whether cached counts need recomputation.
- `requires_prop`: Whether the row still needs propagation processing.
"""
mutable struct XorConstraintMeta
    nnz_par::Int 
    nnz_conj::Int
    is_pure_xor::Bool
    requires_update::Bool
    requires_prop::Bool
end

"""
    XorConstraint

Represent one parity row in XOR or XOR-AND form.

The row stores parity literals in `par`, optional pairwise conjunction terms in
`conj`, and a Boolean right-hand side. Cached structural metadata lives in
`meta`.

# Fields
- `par`: Active parity support as a bit vector.
- `conj`: Symmetric conjunction matrix, or `nothing` for pure XOR rows.
- `rhs`: Right-hand side parity bit.
- `meta`: Cached row statistics and propagation flags.
"""
mutable struct XorConstraint
    par::BitVector
    conj::Union{BitMatrix, Nothing}
    rhs::Bool
    meta::XorConstraintMeta
end

@inline function _count_conj_terms(conj::BitMatrix)
    conj_entries = count(conj)
    @assert iseven(conj_entries)
    return conj_entries ÷ 2
end

@inline _xor_meta(nnz_par::Int, nnz_conj::Int) = XorConstraintMeta(nnz_par, nnz_conj, nnz_conj == 0, false, true)

"""
    XorConstraint(par, conj, rhs)

Construct an XOR-AND constraint from parity and conjunction supports.

# Arguments
- `par`: Parity support bit vector.
- `conj`: Symmetric conjunction matrix with the same dimension as `par`.
- `rhs`: Right-hand side parity bit.

# Returns
- An `XorConstraint` with initialized metadata.

# Throws
- `AssertionError` if `conj` does not match the size of `par`.
"""
function XorConstraint(par::BitVector, conj::BitMatrix, rhs::Bool)
    n = length(par)
    @assert size(conj) == (n, n)
    nnz_conj = _count_conj_terms(conj)
    return XorConstraint(par, nnz_conj == 0 ? nothing : conj, rhs, _xor_meta(count(par), nnz_conj))
end

"""
    XorConstraint(par, rhs)

Construct a pure XOR constraint.
"""
XorConstraint(par::BitVector, rhs::Bool) = XorConstraint(par, nothing, rhs, _xor_meta(count(par), 0))

@inline function _mark_dirty!(con::XorConstraint)
    con.meta.requires_update = true
    con.meta.requires_prop = true
    return con
end


"""
    update!(con)

Refresh cached metadata for `con`.

Recompute the parity and conjunction counts, update the pure-XOR flag, and drop
the conjunction matrix when it becomes empty.

# Returns
- The mutated `con`.

# Side Effects
- Mutates `con.meta` and may set `con.conj = nothing`.
"""
function update!(con::XorConstraint)
    con.meta.nnz_par = count(con.par)
    con.meta.nnz_conj = con.conj === nothing ? 0 : _count_conj_terms(con.conj)
    con.meta.is_pure_xor = con.meta.nnz_conj == 0
    con.meta.is_pure_xor && (con.conj = nothing)
    con.meta.requires_update = false

    return con
end


"""
    xor_con!(con, other)

XOR `other` into `con` in place.

# Arguments
- `con`: Destination row mutated in place.
- `other`: Source row to XOR into `con`.

# Returns
- The mutated `con`.

# Side Effects
- Mutates the parity support, conjunction support, and right-hand side of `con`.

# Notes
- This assumes `con` and `other` are structurally compatible rows of the same
  dimension.
"""
function xor_con!(con::XorConstraint, other::XorConstraint)
    _mark_dirty!(con)
    con.par .⊻= other.par
    @assert !(con.meta.is_pure_xor && !other.meta.is_pure_xor)
    con.conj !== nothing && other.conj !== nothing && !con.meta.is_pure_xor && !other.meta.is_pure_xor &&
        ((con.conj::BitMatrix) .⊻= (other.conj::BitMatrix))
    con.rhs ⊻= other.rhs

    return con
end


abstract type _TripartitePropagationAction end

"""
    _TripartiteImplicationAction

Store the implication-only propagation derived from a tripartite XOR-AND row.
"""
struct _TripartiteImplicationAction <: _TripartitePropagationAction
    x_idx::Int
    y_idx::Int
    r1::Bool
    r2::Bool
end

"""
    _TripartiteRewriteAction

Store the two XOR rows derived from a tripartite XOR-AND rewrite.
"""
struct _TripartiteRewriteAction <: _TripartitePropagationAction
    support1::BitVector
    rhs1::Bool
    support2::BitVector
    rhs2::Bool
end

@inline function _first_nonzero_conj_column(conj::BitMatrix)
    @inbounds for col in 1:size(conj, 2)
        findfirst(@view conj[:, col]) !== nothing && return col
    end

    return nothing
end

@inline function _first_true_in_column(conj::BitMatrix, col::Int)
    return findfirst(@view conj[:, col])
end

@inline function _column_matches_tripartite_pattern(
    conj::BitMatrix,
    col::Int,
    ref_i::Int,
    ref_j::Int,
    mode::UInt8,
)
    @inbounds for row in 1:size(conj, 1)
        expected = if mode == 0x00
            false
        elseif mode == 0x01
            conj[row, ref_i]
        elseif mode == 0x02
            conj[row, ref_j]
        else
            conj[row, ref_i] ⊻ conj[row, ref_j]
        end
        conj[row, col] == expected || return false
    end

    return true
end

"""
    _tripartite_action(con)

Detect a tripartite propagation pattern in `con`.

Inspect an up-to-date XOR-AND row and identify whether it matches the
tripartite implication or rewrite patterns used by parity propagation.

# Arguments
- `con`: Up-to-date XOR-AND constraint to analyze.

# Returns
- `nothing` when no valid tripartite action applies.
- A `_TripartiteImplicationAction` for the implication-only case.
- A `_TripartiteRewriteAction` containing the two derived XOR rows.

# Throws
- `AssertionError` if `con` is not up to date or if multiple valid tripartite
  interpretations are found.
"""
function _tripartite_action(con::XorConstraint)
    @assert !con.meta.requires_update
    con.meta.nnz_conj == 0 && return nothing
    con.conj === nothing && return nothing

    i = _first_nonzero_conj_column(con.conj::BitMatrix)
    i === nothing && return nothing
    j = _first_true_in_column(con.conj::BitMatrix, i)
    @assert j !== nothing

    n = length(con.par)
    nx = 0
    ny = 0
    nz = 0
    x_idx = 0
    y_idx = 0

    @inbounds for v in 1:n
        in_i = (con.conj::BitMatrix)[v, i]
        in_j = (con.conj::BitMatrix)[v, j]

        if in_i
            if in_j
                nz += 1
                _column_matches_tripartite_pattern(con.conj::BitMatrix, v, i, j, 0x03) || return nothing
            else
                ny += 1
                y_idx = v
                _column_matches_tripartite_pattern(con.conj::BitMatrix, v, i, j, 0x02) || return nothing
            end
        elseif in_j
            nx += 1
            x_idx = v
            _column_matches_tripartite_pattern(con.conj::BitMatrix, v, i, j, 0x01) || return nothing
        else
            _column_matches_tripartite_pattern(con.conj::BitMatrix, v, i, j, 0x00) || return nothing
        end
    end

    found_pair = false
    valid_r1 = false
    valid_r2 = false

    for r1 in (false, true), r2 in (false, true)
        matches = true

        @inbounds for v in 1:n
            expected = ((con.conj::BitMatrix)[v, i] & (con.conj::BitMatrix)[v, j]) ⊻
                (r2 & (con.conj::BitMatrix)[v, j]) ⊻
                (r1 & (con.conj::BitMatrix)[v, i])
            con.par[v] == expected || begin
                matches = false
                break
            end
        end

        matches || continue
        ((r1 & r2) == con.rhs && (nx > 1 || ny > 1 || nz > 0)) && continue

        @assert !found_pair "multiple valid (r1, r2) combinations found"
        found_pair = true
        valid_r1 = r1
        valid_r2 = r2
    end

    found_pair || return nothing

    if (valid_r1 & valid_r2) == con.rhs
        @assert nx == 1 && ny == 1 && nz == 0
        @assert x_idx != 0
        @assert y_idx != 0
        return _TripartiteImplicationAction(x_idx, y_idx, valid_r1, valid_r2)
    end

    return _TripartiteRewriteAction(
        copy(@view (con.conj::BitMatrix)[:, j]),
        !valid_r1,
        copy(@view (con.conj::BitMatrix)[:, i]),
        !valid_r2,
    )
end

"""
    fix_var!(con, var_idx, val)

Fix `x[var_idx]` to `val` in place.

# Arguments
- `con`: Constraint mutated in place.
- `var_idx`: Variable position in row coordinates.
- `val`: Fixed Boolean value.

# Returns
- The mutated `con`.

# Side Effects
- Removes all occurrences of `var_idx` from the row and updates the right-hand
  side as needed.
"""
function fix_var!(con::XorConstraint, var_idx::Int, val::Bool)
    @assert 1 <= var_idx <= length(con.par)

    if con.par[var_idx]
        con.rhs ⊻= val
        con.par[var_idx] = false
    end

    if con.conj !== nothing
        conj = con.conj::BitMatrix
        @inbounds for k in eachindex(con.par)
            conj[k, var_idx] || continue
            val && (con.par[k] ⊻= true)
            conj[k, var_idx] = false
            conj[var_idx, k] = false
        end
    end

    return _mark_dirty!(con)
end

"""
    _rewrite_conjunctive_neighbors!(con, conj, var_idx, subst_mask, neg)

Rewrite conjunction terms incident to `var_idx`.

Apply the conjunctive part of a substitution
`x[var_idx] <- (⊻_{j : subst_mask[j]} x[j]) ⊻ neg` while preserving symmetry of
`conj`.

# Returns
- `true` if any conjunctive neighbor was rewritten.
- `false` if `var_idx` had no conjunctive neighbors.

# Side Effects
- Mutates `con.par` and `conj`.
"""
function _rewrite_conjunctive_neighbors!(
    con::XorConstraint,
    conj::BitMatrix,
    var_idx::Int,
    subst_mask::BitVector,
    neg::Bool,
)
    changed = false

    @inbounds for k in eachindex(con.par)
        conj[k, var_idx] || continue
        changed = true
        neg && (con.par[k] ⊻= true)

        for j in eachindex(subst_mask)
            subst_mask[j] || continue
            if j == k
                con.par[j] ⊻= true
            else
                conj[j, k] ⊻= true
                conj[k, j] = conj[j, k]
            end
        end

        conj[k, var_idx] = false
        conj[var_idx, k] = false
    end

    return changed
end

"""
    substitute_var!(con, var_idx, subst_idx, neg)

Substitute one row variable by another row variable, optionally negated.

Apply `x[var_idx] <- x[subst_idx] ⊻ neg` to both parity and conjunction terms of
`con`.

# Arguments
- `con`: Constraint mutated in place.
- `var_idx`: Variable position to eliminate.
- `subst_idx`: Replacement variable position.
- `neg`: Whether to negate the replacement literal.

# Returns
- The mutated `con`.

# Side Effects
- Mutates `con.par`, `con.conj`, and `con.rhs`.

# Notes
- Assumes `con.conj` is symmetric with a false diagonal when present.

# Throws
- `AssertionError` if the supplied indices are invalid or identical.
"""
function substitute_var!(con::XorConstraint, var_idx::Int, subst_idx::Int, neg::Bool)
    @assert 1 <= var_idx <= length(con.par)
    @assert 1 <= subst_idx <= length(con.par)
    @assert var_idx != subst_idx

    @assert con.conj === nothing || con.conj[var_idx, var_idx] == false
    @assert con.conj === nothing || con.conj[subst_idx, subst_idx] == false

    old_par = con.par[var_idx]

    # linear term: old_par * x_var -> old_par * (x_subst ⊻ neg)
    con.par[subst_idx] ⊻= old_par
    con.rhs ⊻= (old_par & neg)
    con.par[var_idx] = false

    con.conj === nothing && return _mark_dirty!(con)

    conj = con.conj::BitMatrix
    @inbounds for k in eachindex(con.par)
        conj[k, var_idx] || continue

        # quadratic terms involving var_idx:
        # (x_var ∧ x_k) -> (x_subst ∧ x_k) ⊻ (neg ∧ x_k)
        conj[k, subst_idx] ⊻= true
        conj[subst_idx, k] = conj[k, subst_idx]
        neg && (con.par[k] ⊻= true)
        conj[k, var_idx] = false
        conj[var_idx, k] = false
    end

    # reduce diagonal: x_j ∧ x_j = x_j
    if conj[subst_idx, subst_idx]
        con.par[subst_idx] ⊻= true
        conj[subst_idx, subst_idx] = false
    end

    return _mark_dirty!(con)
end


"""
    substitute_var!(con, var_idx, subst_mask, neg)

Substitute one row variable by an XOR of several row variables.

Apply

    x[var_idx] <- (⊻_{j : subst_mask[j]} x[j]) ⊻ neg

to both parity and conjunction terms of `con`.

# Arguments
- `con`: Constraint mutated in place.
- `var_idx`: Variable position to eliminate.
- `subst_mask`: Boolean mask selecting replacement variables.
- `neg`: Whether to negate the replacement XOR.

# Returns
- The mutated `con`.

# Side Effects
- Mutates `con.par`, `con.conj`, and `con.rhs`.

# Notes
- Assumes `con.conj` is symmetric with a false diagonal when present.
- `subst_mask[var_idx]` must be `false`.
"""
function substitute_var!(
    con::XorConstraint,
    var_idx::Int,
    subst_mask::BitVector,
    neg::Bool,
)
    @assert 1 <= var_idx <= length(con.par)
    @assert length(subst_mask) == length(con.par)
    @assert !subst_mask[var_idx]

    old_par = con.par[var_idx]

    # linear term:
    # old_par * x_i -> old_par * ((xor-sum over subst vars) ⊻ neg)
    if old_par
        con.par .⊻= subst_mask
        con.rhs ⊻= neg
    end
    con.par[var_idx] = false

    con.conj === nothing && return _mark_dirty!(con)

    conj = con.conj::BitMatrix
    _rewrite_conjunctive_neighbors!(con, conj, var_idx, subst_mask, neg)
    return _mark_dirty!(con)
end

"""
    substitute_var_in_conjunctive_terms!(con, var_idx, subst_mask, neg)

Rewrite only the conjunction terms incident to `var_idx`.

Apply

    x[var_idx] <- (⊻_{j : subst_mask[j]} x[j]) ⊻ neg

inside conjunction terms of `con`, while leaving any linear occurrence of
`x[var_idx]` untouched.

# Arguments
- `con`: Constraint mutated in place.
- `var_idx`: Variable position to rewrite in conjunction terms.
- `subst_mask`: Boolean mask selecting replacement variables.
- `neg`: Whether to negate the replacement XOR.

# Returns
- The mutated `con`.

# Side Effects
- Mutates `con.conj` and possibly `con.par`.

# Notes
- Assumes `con.conj` is symmetric with a false diagonal when present.
- `subst_mask[var_idx]` must be `false`.
"""
function substitute_var_in_conjunctive_terms!(
    con::XorConstraint,
    var_idx::Int,
    subst_mask::BitVector,
    neg::Bool,
)
    @assert 1 <= var_idx <= length(con.par)
    @assert length(subst_mask) == length(con.par)
    @assert !subst_mask[var_idx]

    con.conj === nothing && return con

    conj = con.conj::BitMatrix
    _rewrite_conjunctive_neighbors!(con, conj, var_idx, subst_mask, neg) || return con
    return _mark_dirty!(con)
end

"""
    XorConstraintBuilder

Build an XOR or XOR-AND constraint in variable-id space.

Terms are toggled modulo 2 in the internal dictionaries and materialized into
bit-vectors and matrices by [`build`](@ref).

# Fields
- `par`: Active parity terms keyed by variable id.
- `conj`: Active conjunction terms keyed by ordered variable-id pairs.
- `rhs`: Right-hand side parity bit.
"""
mutable struct XorConstraintBuilder
    par::Dict{VarId, Bool}
    conj::Dict{Tuple{VarId, VarId}, Bool}
    rhs::Bool
end

XorConstraintBuilder() = XorConstraintBuilder(Dict{VarId, Bool}(), Dict{Tuple{VarId, VarId}, Bool}(), false)

"""
    add_par!(builder, var_id)

Toggle linear term `var_id` in `builder`.

# Returns
- The new active/inactive state stored for `var_id`.
"""
add_par!(builder::XorConstraintBuilder, var_id::VarId) = (builder.par[var_id] = !get(builder.par, var_id, false))

"""
    add_conj!(builder, var_id_1, var_id_2)

Toggle conjunctive term `(var_id_1, var_id_2)` in `builder`.

# Returns
- The new active/inactive state stored for the ordered variable pair.
"""
function add_conj!(builder::XorConstraintBuilder, var_id_1::VarId, var_id_2::VarId)
    var_id_1 > var_id_2 && ((var_id_1, var_id_2) = (var_id_2, var_id_1))
    return builder.conj[(var_id_1, var_id_2)] = !get(builder.conj, (var_id_1, var_id_2), false)
end

"""
    negate!(builder)

Toggle the builder right-hand side bit.

# Returns
- The new right-hand side bit.
"""
negate!(builder::XorConstraintBuilder) = (builder.rhs ⊻= true)

"""
    build(builder, nvars, var_to_pos)

Materialize `builder` into an `XorConstraint` using `var_to_pos`.

# Arguments
- `builder`: Constraint builder in variable-id space.
- `nvars`: Number of parity variables in the target row space.
- `var_to_pos`: Mapping from variable ids to row positions.

# Returns
- An `XorConstraint` over row positions `1:nvars`.
"""
function build(builder::XorConstraintBuilder, nvars::Int, var_to_pos::Dict{VarId, Int})
    par = falses(nvars)
    conj = falses(nvars, nvars)
    is_pure_xor = true

    for (var_id, val) in builder.par
        val && (par[var_to_pos[var_id]] = true)
    end

    for ((var_id_1, var_id_2), val) in builder.conj
        conj[var_to_pos[var_id_1], var_to_pos[var_id_2]] = val
        conj[var_to_pos[var_id_2], var_to_pos[var_id_1]] = val
        is_pure_xor = false
    end

    con = if is_pure_xor
        XorConstraint(par, builder.rhs)
    else
        XorConstraint(par, conj, builder.rhs)
    end

    return con
end
