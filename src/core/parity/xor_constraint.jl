const VarId = Int


"""
struct for storing XOR constraints
"""
mutable struct XorConstraint
    pos_to_var::Vector{VarId}
    var_to_pos::Dict{VarId, Int}
    par::BitVector
    conj::BitMatrix
    rhs::Bool
end

"""
xor other onto con 
"""
function xor_con!(con::XorConstraint, other::XorConstraint)
    con.par .⊻= other.par
    con.conj .⊻= other.conj
    con.rhs ⊻= other.rhs
end




"""
Builder struct for constructing XOR-AND constraints

"""
mutable struct XorConstraintBuilder
    par::Dict{VarId, Bool}
    conj::Dict{Tuple{VarId, VarId}, Bool}
    rhs::Bool
end

XorConstraintBuilder() = XorConstraintBuilder(Dict{VarId, Bool}(), Dict{Tuple{VarId, VarId}, Bool}(), false) 

function add_par!(builder::XorConstraintBuilder, var_id::VarId)
    builder.par[var_id] = !get(builder.par, var_id, false)
end

function add_conj!(builder::XorConstraintBuilder, var_id_1::VarId, var_id_2::VarId)
    if var_id_1 > var_id_2 
        var_id_1, var_id_2 = var_id_2, var_id_1
    end 
    builder.conj[(var_id_1, var_id_2)] = !get(builder.conj, (var_id_1, var_id_2), false)
end

negate!(builder::XorConstraintBuilder) = (builder.rhs ⊻= true)

function build(builder::XorConstraintBuilder, nvars::Int, pos_to_var::Vector{VarId}, var_to_pos::Dict{VarId, Int})
    par = falses(nvars)
    conj = falses(nvars, nvars)

    for (var_id, val) in builder.par
        val && (par[var_to_pos[var_id]] = true)
    end

    for ((var_id_1, var_id_2), val) in builder.conj
        conj[var_to_pos[var_id_1], var_to_pos[var_id_2]] = val
        conj[var_to_pos[var_id_2], var_to_pos[var_id_1]] = val
    end

    con = XorConstraint(pos_to_var, var_to_pos, par, conj, builder.rhs)

    return con
end

