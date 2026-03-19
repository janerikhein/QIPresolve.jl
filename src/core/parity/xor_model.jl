

mutable struct ParityModel

    #TODO: store parity and conjunctioon pivots
    var_id_to_pos::Dict{VarId, Int}
    pos_to_var_id::Vector{Int}
    cons::Vector{XorConstraint}
end

function ParityModel(var_to_pos::Dict{VarId, Int}, pos_to_var::Vector{Int}, cons::Vector{XorConstraint})
    #TODO: implement constructor
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