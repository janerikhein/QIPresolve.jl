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
function build(builder::XorConstraintBuilder, nvars::Int, pos_to_var::Vector{VarId}, var_to_pos::Dict{VarId, Int})
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
