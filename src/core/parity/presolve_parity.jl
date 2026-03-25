

function build_parity_model(model::QPModel)

    con_builders = XorConstraintBuilder[]

    for con in model.cons
        
        con_vars = vars(con) 

        # Build Mod-2 XOR constraint and check if mod-4 relaxation is applicable
        builder_mod2 = XorConstraintBuilder()
        discard_mod_4 = false

        for vid in con_vars
            is_bin = is_binary(model.vars[vid])
            diag = convert(Int, get_quad_coeff(con.qe, vid, vid))
            lin = convert(Int, get_lin_coeff(con.qe, vid))
            # mod 4: check if mod 4 is applicable, i.e. no non-binary odd linear coefficients
            if !is_bin && mod(lin, 2) == 1
                discard_mod_4 = true
            end
            # mod 2: check for odd (diagonal + linear)
            if mod(diag + lin, 2) == 1
                add_par!(builder_mod2, vid)
            end
        end

        push!(con_builders, builder_mod2)

        # Build mod-4 XOR-AND constraint if applicable
        discard_mod_4 && continue
        builder_mod4 = XorConstraintBuilder()

        for (i, vid_i) in enumerate(con_vars)
            is_bin_i = is_binary(model.vars[vid_i])
            diag_i = convert(Int, get_quad_coeff(con.qe, vid_i, vid_i))
            lin_i = convert(Int, get_lin_coeff(con.qe, vid_i))
            @assert mod(lin_i, 2) == 0

            # mod 4: add parity terms from non-binary linear terms with even coefficients
            if !is_bin_i && mod(lin_i ÷ 2, 2) == 1
                add_par!(builder_mod4, var_id_i)
            end

            # mod 4: add parity term from (diagonal + binary linear) terms equal 2 or 3 (mod 4)
            ci = is_bin_i ? diag_i + lin_i : diag_i
            if mod(ci, 4) >> 1 == 1
                add_par!(builder_mod4, var_id_i)
            end

            for j in i+1:lastindex(con_vars)
                vid_j = con_vars[j]
                is_bin_j = is_binary(model.vars[vid_i])
                diag_j = convert(Int, get_quad_coeff(con.qe, vid_j, vid_j))
                lin_j = convert(Int, get_lin_coeff(con.qe, vid_j))
                bilin_ij = convert(Int, get_quad_coeff(con.qe, vid_i, vid_j))
                #TODO
            end


    end
   
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