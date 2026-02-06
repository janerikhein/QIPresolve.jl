const VarId = Int
const INVALID_VAR_ID = -1

struct IndexRef
    id1::VarId
    id2::VarId
    quad::Bool
end

Base.:(==)(a::IndexRef, b::IndexRef) =
    a.id1 == b.id1 &&
    a.id2 == b.id2 &&
    a.quad == b.quad

Base.hash(r::IndexRef, h::UInt) =
    hash(r.id1, hash(r.id2, hash(r.quad, h)))

function Base.show(io::IO, r::IndexRef)
    if r.quad
        print(io, "p_$(r.id1)∧p_$(r.id2)")
    else
        print(io, "p_$(r.id1)")
    end
end

QuadRef(id1::VarId, id2::VarId) = IndexRef(id1, id2, true)
LinRef(id::VarId) = IndexRef(id1, INVALID_VAR_ID, false)


mutable struct XorAndModel
    xor_matrix::Vector{BitVector}
    idx_to_var_id::Vector{IndexRef}
    var_id_to_idx::Dict{IndexRef, Int}
end

struct ParityState
    model::QPModel
    xor_model::XorAndModel
end

function presolve!(presolve_state::ParityState)
    xor_mat = presolve_state.xor_model.xor_matrix
    gauss_jordan_gf2!(xor_mat, FirstOnePivot())
    for vec in xor_mat
        for i in findall(vec)
            
        end
    end

end