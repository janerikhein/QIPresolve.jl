const VarId = Int


"""
    - var_id_to_pos
    - pos_to_var_id
    - cons::Vector{XorConstraint} 
        main matrix representation + rhs
    - propagator::ImplicationNetwork.ImpGraph
        ref to implication network for propagation
    - pivots::Vector{Int}
        pivots[i] = j, if (i,j) is pivot element, or pivots[i] = 0 if i not pivot row
"""
mutable struct ParityModel
    var_id_to_pos::Dict{VarId, Int}
    pos_to_var_id::Vector{Int}
    cons::Vector{XorConstraint}
    propagator::ImpGraph
end

function ParityModel(var_to_pos::Dict{VarId, Int}, pos_to_var::Vector{Int}, cons::Vector{XorConstraint})
    nvars = length(pos_to_var)
    implication_graph = ImpGraph(nvars)
    return ParityModel(var_to_pos, pos_to_var, cons, implication_graph)
end

function get_parity_model(model::QPModel)

    # build relevant variables ids and constraint indices
    con_indices = Int[]
    var_id_set = Set{Int}()
    for (i, con) in enumerate(model.cons)
        # skip inequalities and constraints with fractional coefficients
        (!is_equality(con) || !is_integer(con)) && continue
        push!(con_indices, i)
        for var_id in vars(con)
            push!(var_id_set, var_id)
        end
    end
    var_ids = collect(var_id_set)

    # create xor constraint builders and count variable occurences in xor cons
    var_counts = Dict{}
    con_builders = XorConstraintBuilder[]
    for idx in con_indices
        con = model.cons[idx]
        builder_mod2 = XorConstraintBuilder()
        builder_mod4 = XorConstraintBuilder()
        discard_mod_4 = false

        for (i, var_id_i) in enumerate(var_ids)
            is_bin_i = is_binary(model.vars[var_id_i])
            diag_i = convert(Int, get_quad_coeff(con.qe, var_id_i, var_id_i))
            lin_i = convert(Int, get_lin_coeff(con.qe, var_id_i))

            # mod 2: check for odd (diagonal + linear)
            if mod(diag_i + lin_i, 2) == 1
                add_par!(builder_mod2, var_id_i)
            end

            # mod 4: check if mod 4 is applicable, i.e. no non-binary odd linear coefficients
            if !is_bin_i && mod(lin_i, 2) == 1
                discard_mod_4 = true
            end
            discard_mod_4 && continue

            # mod 4: add parity terms from non-binary linear terms with even coefficients
            if !is_bin_i && mod(lin_i, 2) == 0 && mod(lin_i >> 1, 2) == 1
                add_par!(builder_mod4, var_id_i)
            end

            # mod 4: add parity term from (diagonal + binary linear) terms equal 2 or 3 (mod 4)
            ci = is_bin_i ? diag_i + lin_i : diag_i
            if mod(ci, 4) >> 1 == 1
                add_par!(builder_mod4, var_id_i)
            end

            for j in (i + 1):lastindex(var_ids)
                var_id_j = var_ids[j]
                is_bin_j = is_binary(model.vars[var_id_j])
                diag_j = convert(Int, get_quad_coeff(con.qe, var_id_j, var_id_j))
                lin_j = convert(Int, get_lin_coeff(con.qe, var_id_j))
                bilin_ij = convert(Int, get_quad_coeff(con.qe, var_id_i, var_id_j))

                # mod 4: check if bilinear term 2qij has odd qij
                if mod(bilin_ij, 4) >> 1 == 1
                    add_conj!(builder_mod4, var_id_i, var_id_j)
                end

                # mod 4: check if both (diagonal + binary linear) equal 1 or 3 (mod 4)
                cj = is_bin_j ? diag_j + lin_j : diag_j
                if mod(ci, 2) == 1 && mod(cj, 2) == 1
                    add_conj!(builder_mod4, var_id_i, var_id_j)
                end
            end
        end

        # set rhs
        con_rhs = convert(Int, con.rhs)
        con_rhs % 2 == 1 && negate!(builder_mod2)
        (con_rhs % 4) >> 1 == 1 && negate!(builder_mod4)

        push!(con_builders, builder_mod2)
        !discard_mod_4 && push!(con_builders, builder_mod4)
    end

    # compute variable ordering based on decreasing variable var count
    var_counts = Dict(vid => 0 for vid in var_ids)
    for builder in con_builders
        for (vid, val) in builder.par
            val && (var_counts[vid] += 1)
        end
        for ((vid1, vid2), val) in builder.conj
            if val
                var_counts[vid1] += 1
                var_counts[vid2] += 1
            end
        end
    end

    # filter out vars with zero count
    pos_to_var = sort!([(vid, c) for (vid, c) in var_counts if c > 0], by = last, rev = true) .|> first
    var_to_pos = Dict(vid => i for (i, vid) in enumerate(pos_to_var))
    nvars = length(pos_to_var)

    # build parity model
    cons = [build(builder, nvars, pos_to_var, var_to_pos) for builder in con_builders]

    return ParityModel(var_to_pos, pos_to_var, cons)
end


function propagate!(model::ParityModel)
    for (i, con) in enumerate(model.cons)
        bipartite_split = split_bipartite(con)
        if bipartite_split !== nothing
            con1, con2 = bipartite_split
            model.cons[i] = con1
            push!(model.cons, con2)
        end
    end
end

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
"""
function apply_substitution!(model::XorAndModel, xor_row_idx::Int, pivot_col_idx::Int)

end
"""

