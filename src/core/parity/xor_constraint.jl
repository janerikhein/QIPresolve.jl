mutable struct XorConstraintMeta
    nnz_par::Int #non-zeros in parity vector
    nnz_conj::Int #non zeros in conjunction matrix
    is_pure_xor::Bool
    requires_update::Bool
    requires_prop::Bool
end

mutable struct XorConstraint
    par::BitVector
    conj::Union{BitMatrix, Nothing}
    rhs::Bool
    meta::XorConstraintMeta
end

function XorConstraint(par::BitVector, conj::BitMatrix, rhs::Bool)
    n = length(par)
    @assert size(conj, 1) == n
    @assert size(conj, 2) == n

    nnz_par = count(par)
    conj_entries = count(conj)
    @assert iseven(conj_entries)
    nnz_conj = conj_entries ÷ 2
    is_pure_xor = nnz_conj == 0
    normalized_conj = is_pure_xor ? nothing : conj
    meta = XorConstraintMeta(nnz_par, nnz_conj, is_pure_xor, false, true)

    return XorConstraint(par, normalized_conj, rhs, meta)
end

function XorConstraint(par::BitVector, rhs::Bool)
    meta = XorConstraintMeta(count(par), 0, true, false, true)
    return XorConstraint(par, nothing, rhs, meta)
end


function update!(con::XorConstraint)
    con.meta.nnz_par = count(con.par)
    con.meta.nnz_conj = 0
    con.meta.is_pure_xor = true

    if con.conj !== nothing
        conj_entries = count(con.conj)
        @assert iseven(conj_entries)
        con.meta.nnz_conj = conj_entries ÷ 2

        if con.meta.nnz_conj == 0
            con.conj = nothing
        else
            con.meta.is_pure_xor = false
        end
    end

    con.meta.requires_update = false

    return con
end


"""
    xor_con!(con::XorConstraint, other::XorConstraint) -> XorConstraint

XOR `other` into `con` in place.
"""
function xor_con!(con::XorConstraint, other::XorConstraint)
    con.meta.requires_update = true
    con.meta.requires_prop = true
    con.par .⊻= other.par
    @assert !(con.meta.is_pure_xor && !other.meta.is_pure_xor)
    if con.conj !== nothing && other.conj !== nothing && !con.meta.is_pure_xor && !other.meta.is_pure_xor
        con.conj .⊻= other.conj
    end
    con.rhs ⊻= other.rhs

    return con
end


"""
    split_bipartite(con::XorConstraint) -> Union{Nothing,Tuple{XorConstraint,XorConstraint}}

Detect whether `con` encodes a full bipartite conjunctive structure and split
it into two pure XOR constraints.

Returns `nothing` when splitting is not applicable.

TODO: implement extended bipartite check
    1) check full bipartite of {(i,j): p_ij = 1, p_i=0, p_j=0} -> get two parts A,B
    2) Let P = {i: p_i=1}
    3) Check if p_ix = 1, iff i in P and x in A∪B
    4) build xor cons over P∪A and P∪B
"""
function split_bipartite(con::XorConstraint)

    # rhs must be 1 for splitting
    !con.rhs && return nothing

    # constraint must be purely conjunctive
    (con.conj === nothing || con.meta.is_pure_xor) && return nothing
    any(con.par) && return nothing

    # degree scan: allow 0 (isolated) plus at most two positive degrees
    n = size(con.conj, 2)
    degree_1 = 0
    degree_2 = 0
    piv = 0
    degrees = Vector{Int}(undef, n)

    @inbounds for v in 1:n
        deg = degrees[v] = count(@view con.conj[:, v])
        if deg == 0
            continue
        elseif degree_1 == 0
            degree_1 = deg
            piv = v
        elseif deg == degree_1
            continue
        elseif degree_2 == 0
            degree_2 = deg
        elseif deg == degree_2
            continue
        else
            return nothing
        end
    end

    # constraint is empty
    piv == 0 && return nothing

    # construct the two possible neighborhoods
    N1 = @view con.conj[:, piv]
    first_neighbor = findfirst(N1)
    N2 = @view con.conj[:, first_neighbor]

    # check if neighborhoods are disjoint
    any(N1 .& N2) && return nothing

    # check for full bipartite
    @inbounds for v in 1:n
        col = @view con.conj[:, v]

        if degrees[v] == 0
            continue  # isolated is allowed
        end

        if N1[v]
            # v is in N1-side => its neighborhood must be exactly N2
            any(col .⊻ N2) && return nothing
        elseif N2[v]
            # v is in N2-side => its neighborhood must be exactly N1
            any(col .⊻ N1) && return nothing
        else
            # v has edges but is in neither N1 nor N2 => impossible
            return nothing
        end
    end

    # build the two resulting XOR constraints from the two neighborhoods
    con1 = XorConstraint(BitVector(copy(N1)), true)
    con2 = XorConstraint(BitVector(copy(N2)), true)

    return con1, con2
end

"""
    fix_var!(con::XorConstraint, var_idx::Int, val::Bool) -> XorConstraint

Fix `x[var_idx]` to `val` in place.
"""
function fix_var!(con::XorConstraint, var_idx::Int, val::Bool)
    @assert 1 <= var_idx <= length(con.par)

    if con.par[var_idx]
        con.rhs ⊻= val
        con.par[var_idx] = false
    end

    if con.conj !== nothing
        conj = con.conj::BitMatrix
        conj_var = BitVector(copy(@view conj[:, var_idx]))

        if val
            con.par .⊻= conj_var
        end

        conj[:, var_idx] .= false
        conj[var_idx, :] .= false
    end

    con.meta.requires_update = true
    con.meta.requires_prop = true

    return con
end

"""
Performs substitution x[var_idx] <- x[subst_idx] ⊻ neg.

Assumes:
- con.conj is symmetric
- con.conj[i, i] == false on entry
- var_idx != subst_idx
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

    if con.conj === nothing
        con.meta.requires_update = true
        con.meta.requires_prop = true
        return con
    end

    conj = con.conj::BitMatrix
    conj_var = view(conj, :, var_idx)
    conj_subst = view(conj, :, subst_idx)

    # quadratic terms involving var_idx:
    # (x_var ∧ x_k) -> (x_subst ∧ x_k) ⊻ (neg ∧ x_k)
    conj_subst .⊻= conj_var

    if neg
        con.par .⊻= conj_var
    end

    # restore symmetry from updated subst column
    conj[subst_idx, :] .= conj_subst

    # remove var_idx completely
    conj_var .= false
    conj[var_idx, :] .= false

    # reduce diagonal: x_j ∧ x_j = x_j
    if conj[subst_idx, subst_idx]
        con.par[subst_idx] ⊻= true
        conj[subst_idx, subst_idx] = false
    end

    con.meta.requires_update = true
    con.meta.requires_prop = true

    return con
end


"""
Performs substitution

    x[var_idx] <- (⊻_{j : subst_mask[j]} x[j]) ⊻ neg

Assumes:
- con.conj is symmetric
- con.conj[i, i] == false on entry for all i
- subst_mask[var_idx] == false
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

    if con.conj === nothing
        con.meta.requires_update = true
        con.meta.requires_prop = true
        return con
    end

    conj = con.conj::BitMatrix
    subst_indices = findall(subst_mask)
    neighbor_indices = findall(@view conj[:, var_idx])

    # quadratic terms:
    # (x_i ∧ x_k) -> ((xor-sum over subst vars) ∧ x_k) ⊻ (neg ∧ x_k)
    for k in neighbor_indices
        if neg
            con.par[k] ⊻= true
        end

        for j in subst_indices
            if j == k
                con.par[j] ⊻= true
            else
                conj[j, k] ⊻= true
                conj[k, j] = conj[j, k]
            end
        end
    end

    # remove var_idx completely
    conj[:, var_idx] .= false
    conj[var_idx, :] .= false

    con.meta.requires_update = true
    con.meta.requires_prop = true

    return con
end

"""
Performs substitution inside conjunction terms only

    x[var_idx] <- (⊻_{j : subst_mask[j]} x[j]) ⊻ neg

This only rewrites terms of the form `x[var_idx] ∧ x[k]` and leaves any linear
occurrence of `x[var_idx]` untouched.

Assumes:
- con.conj is symmetric
- con.conj[i, i] == false on entry for all i
- subst_mask[var_idx] == false
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
    subst_indices = findall(subst_mask)
    neighbor_indices = findall(@view conj[:, var_idx])
    isempty(neighbor_indices) && return con

    # quadratic terms:
    # (x_i ∧ x_k) -> ((xor-sum over subst vars) ∧ x_k) ⊻ (neg ∧ x_k)
    for k in neighbor_indices
        if neg
            con.par[k] ⊻= true
        end

        for j in subst_indices
            if j == k
                con.par[j] ⊻= true
            else
                conj[j, k] ⊻= true
                conj[k, j] = conj[j, k]
            end
        end
    end

    # remove var_idx from conjunction terms
    conj[:, var_idx] .= false
    conj[var_idx, :] .= false

    con.meta.requires_update = true
    con.meta.requires_prop = true

    return con
end

"""
    XorConstraintBuilder

Builder for XOR/XOR-AND constraints in variable-id space.

Terms are toggled modulo 2 in the internal dictionaries and materialized into
bit-vectors/matrices by [`build`](@ref).
"""
mutable struct XorConstraintBuilder
    par::Dict{VarId, Bool}
    conj::Dict{Tuple{VarId, VarId}, Bool}
    rhs::Bool
end

"""
    XorConstraintBuilder()

Create an empty constraint builder with `rhs == false`.
"""
XorConstraintBuilder() = XorConstraintBuilder(Dict{VarId, Bool}(), Dict{Tuple{VarId, VarId}, Bool}(), false)

"""
    add_par!(builder::XorConstraintBuilder, var_id::VarId) -> Bool

Toggle linear term `var_id` in `builder`.
"""
function add_par!(builder::XorConstraintBuilder, var_id::VarId)
    return builder.par[var_id] = !get(builder.par, var_id, false)
end

"""
    add_conj!(builder::XorConstraintBuilder, var_id_1::VarId, var_id_2::VarId) -> Bool

Toggle conjunctive term `(var_id_1, var_id_2)` in `builder`.
"""
function add_conj!(builder::XorConstraintBuilder, var_id_1::VarId, var_id_2::VarId)
    if var_id_1 > var_id_2
        var_id_1, var_id_2 = var_id_2, var_id_1
    end
    return builder.conj[(var_id_1, var_id_2)] = !get(builder.conj, (var_id_1, var_id_2), false)
end

"""
    negate!(builder::XorConstraintBuilder) -> Bool

Toggle the builder right-hand side bit.
"""
negate!(builder::XorConstraintBuilder) = (builder.rhs ⊻= true)

"""
    build(builder::XorConstraintBuilder, nvars::Int, pos_to_var, var_to_pos) -> XorConstraint

Materialize `builder` into an `XorConstraint` using the provided variable
position mappings.
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

function build(
    builder::XorConstraintBuilder,
    nvars::Int,
    pos_to_var,
    var_to_pos::Dict{VarId, Int},
)
    return build(builder, nvars, var_to_pos)
end
