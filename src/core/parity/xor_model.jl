

mutable struct ParityModel
    var_id_to_pos::Dict{VarId, Int}
    pos_to_var_id::Vector{Int}
    cons::Vector{XorConstraint}
    propagator::ParityPropagator
end

function ParityModel(var_to_pos::Dict{VarId, Int}, pos_to_var::Vector{Int}, cons::Vector{XorConstraint})
    cons::Vector{XorConstraint}
    propagator::ParityPropagator 
end



function reformulate_bipartite_cons!(model::ParityModel)

end


function has_unpivoted_con(model::ParityModel)

end

function has_unpivoted_xor_con(model::ParityModel)

end


function gauss_jordan_xor!(model::ParityModel)
    
end


function gauss_jordan_xor_and!(model::ParityModel)
    
end


function substitute_parity_pivots!(model::ParityModel)

end
"""
    gauss_jordan!(model::ParityModel; skip_conj=false) -> ParityModel

Perform sparse Gauss-Jordan elimination over `GF(2)` on `model.cons`.

At each step, the sparsest non-redundant row is selected as pivot. If
`skip_conj=true`, rows containing conjunction terms are ignored and elimination
is restricted to pure XOR constraints.
"""
function gauss_jordan!(model::ParityModel; skip_conj::Bool = false)
    pivots = [false for _ in 1:length(model.cons)]
    while true
        piv_row = 0
        min_nnz = typemax(Int)
        for (i, con) in enumerate(model.cons)
            # skip existing pivot rows
            pivots[i] && continue

            # skip xor-and cons if flag is set
            skip_conj && !con.is_pure_xor && continue

            # update con if required
            con.has_changed && update!(con)

            # skip redundant cons
            con.nnz == 0 && continue

            # compute sparsest row and choose as pivot row
            if con.nnz < min_nnz
                piv_row = i
                min_nnz = con.nnz
            end
        end
        # no further pivot row found
        piv_row == 0 && break

        piv_con = model.cons[piv_row]
        pivots[piv_row] = true

        # get pivot element index. always choose from conjunction terms if possible (i.e. not pure xor)
        piv_col_idx1, piv_col_idx2 = 0, 0
        if model.cons[piv_row].is_pure_xor
            piv_col_idx1 = findfirst(piv_con.par)
        else
            @assert skip_conj == false
            @assert any(piv_con.conj)
            piv_col_idx1, piv_col_idx2 = findfirst(piv_con.conj).I
        end

        for (i, con) in enumerate(model.cons)
            i == piv_row && continue

            # skip xor-and cons if flag is set
            skip_conj && !con.is_pure_xor && continue

            # skip redundant cons
            con.nnz == 0 && continue

            # pivot element is parity term
            if piv_col_idx2 == 0
                if con.par[piv_col_idx1]
                    xor_con!(con, piv_con)
                end
            else # pivot element is conjunctive term
                if con.conj !== nothing && con.conj[piv_col_idx1, piv_col_idx2]
                    xor_con!(con, piv_con)
                end
            end
        end
    end

    return model
end