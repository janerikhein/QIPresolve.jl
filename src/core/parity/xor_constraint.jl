const VarId = Int


"""
struct for storing XOR constraints
"""
mutable struct XorConstraint
    pos_to_var::Vector{VarId}
    var_to_pos::Dict{VarId, Int}
    par::BitVector
    conj::Union{BitMatrix, Nothing}
    rhs::Bool
    conj_deg::Union{Vector{Int}, Nothing}
    nnz::Int
    is_pure_xor::Bool
    has_changed::Bool
end

function XorConstraint(pos_to_var::Vector{VarId}, var_to_pos::Dict{VarId, Int}, par::BitVector, conj::BitMatrix, rhs::Bool)
    n = size(conj, 2)
    conj_deg = Vector{Int}(undef, n)
    for i in 1:n
        conj_deg[i] = count(@view conj[:, i])
    end
    conj_deg_sum = sum(conj_deg)
    @assert conj_deg_sum % 2 == 0
    nnz = count(par) + sum(conj_deg_sum) ÷ 2

    return XorConstraint(pos_to_var, var_to_pos, par, conj, rhs, conj_deg, nnz, false, false)
end

function XorConstraint(pos_to_var::Vector{VarId}, var_to_pos::Dict{VarId, Int}, par::BitVector, rhs::Bool)
    nnz = count(par)
    return XorConstraint(pos_to_var, var_to_pos, par, nothing, rhs, nothing, nnz, true, false)
end

function update!(con::XorConstraint)
    conj_deg_sum = 0
    if con.conj !== nothing 
        @inbounds for i in eachindex(con.conj_deg)
            con.conj_deg[i] = count(@view con.conj[:, i])
        end

        for deg in con.conj_deg
            @assert !(deg > 0 && con.is_pure_xor) 
            conj_deg_sum += deg
        end

        @assert conj_deg_sum % 2 == 0
        if conj_deg_sum == 0
            con.is_pure_xor = true
            con.conj = nothing
            con.conj_deg = nothing
        end
    end

    con.nnz = count(con.par) + conj_deg_sum ÷ 2
    con.has_changed = false
end

"""
xor other onto con 
"""
function xor_con!(con::XorConstraint, other::XorConstraint; update::Bool = true)
    con.has_changed = true
    con.par .⊻= other.par
    @assert !(con.is_pure_xor && !other.is_pure_xor)
    if con.conj !== nothing && other.conj !== nothing && !con.is_pure_xor && !other.is_pure_xor
        con.conj .⊻= other.conj
    end
    con.rhs ⊻= other.rhs
    update && update!(con)

    return con
end

"""
    neccessary conditions for constraint to be full bipartite allowing split
"""
function split_bipartite(con::XorConstraint)

    # rhs must be 1 for splitting
    !con.rhs && return nothing

    # constraint must be purely conjunctive
    (con.conj === nothing || con.is_pure_xor) && return nothing
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
    con1 = XorConstraint(con.pos_to_var, con.var_to_pos, copy(N1), true)
    con2 = XorConstraint(con.pos_to_var, con.var_to_pos, copy(N2), true)

    return con1, con2
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
    return builder.par[var_id] = !get(builder.par, var_id, false)
end

function add_conj!(builder::XorConstraintBuilder, var_id_1::VarId, var_id_2::VarId)
    if var_id_1 > var_id_2
        var_id_1, var_id_2 = var_id_2, var_id_1
    end
    return builder.conj[(var_id_1, var_id_2)] = !get(builder.conj, (var_id_1, var_id_2), false)
end

negate!(builder::XorConstraintBuilder) = (builder.rhs ⊻= true)

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
