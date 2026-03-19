function build_parity_model(model::QPModel)

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

                # by assumption bilinear terms are represented as 2qij
                @assert mod(bilin_ij, 2) == 0

                # mod 4: check if bilinear term 2qij has odd qij
                if mod(bilin_ij, 4) != 0
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
        mod(con_rhs, 2) == 1 && negate!(builder_mod2)
        (mod(con_rhs, 4) >> 1) == 1 && negate!(builder_mod4)

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

# TODO: add postsolve stack structure here
function parity_presolve!(model::QPModel)
    
    #TODO: initialize propagator
    propagator = ParityPropagator()

    while true
        # build parity model
        parity_model = build_parity_model(model)
        reformulate_bipartite_cons!(parity_model)
        propagate!(parity_model, propagator)

        while has_unpivoted_con(parity_model)

            # Phase 1: GJ on XOR cons
            if has_unpivoted_xor_con(parity_model)
                gauss_jordan_xor!(parity_model)
                propagate!(parity_model, propagator)
                substitute_parity_pivots!(parity_model)

            # Phase 2: GJ on XOR-AND cons
            else
                gauss_jordan_xor_and!(parity_model)
                propagate!(parity_model, propagator)
            end  
        end

        parities_fixed = fix_parities!(model, propagator) 
        patterns_added = fix_parity_patterns!(model, propagator)
        
        parities_fixed || patterns_added || break
    end

    add_binary_implications!(model, propagator)

end