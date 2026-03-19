

mutable struct XorConstraint
    # index maps
    pos_to_var::Vector{VarId}
    var_to_pos::Dict{VarId, Int}
    # data
    par::BitVector
    conj::Union{BitMatrix, Nothing}
    rhs::Bool
    # meta data
    nnz_par::Int
    nnz_conj::Int
    conj_deg::Union{Vector{Int}, Nothing}
    is_pure_xor::Bool
    last_changed::Int
    last_updated::Int
    last_prop::Int
end

function XorConstraint(pos_to_var::Vector{VarId}, var_to_pos::Dict{VarId, Int}, par::BitVector, conj::BitMatrix, rhs::Bool)
    #TODO: Implement constructor
end


function update!(con::XorConstraint, cur_time::Int)
    
    #TODO: update metadata of con if last_changed > last_updated and set last_updated to cur_time
end


"""
    xor_con!(con::XorConstraint, other::XorConstraint; update=true) -> XorConstraint

XOR `other` into `con` in place.
"""
function xor_con!(con::XorConstraint, other::XorConstraint)
    con.has_changed = true
    con.par .⊻= other.par
    @assert !(con.is_pure_xor && !other.is_pure_xor)
    if con.conj !== nothing && other.conj !== nothing && !con.is_pure_xor && !other.is_pure_xor
        con.conj .⊻= other.conj
    end
    con.rhs ⊻= other.rhs

    return con
end


function split_bipartite(con::XorConstraint)
    #TODO: implement split bipartite for extended case
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
        XorConstraint(pos_to_var, var_to_pos, par, builder.rhs)
    else
        XorConstraint(pos_to_var, var_to_pos, par, conj, builder.rhs)
    end

    return con
end
