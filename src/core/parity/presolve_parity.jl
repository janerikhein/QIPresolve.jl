function get_builders(model::QPModel, con::Constraint)
    con_vars = vars(con)

    # Build mod-2 XOR constraint and check if mod-4 relaxation is applicable.
    builder_mod2 = XorConstraintBuilder()
    discard_mod_4 = false

    for vid in con_vars
        is_bin = is_binary(model.vars[vid])
        diag = convert(Int, get_quad_coeff(con.qe, vid, vid))
        lin = convert(Int, get_lin_coeff(con.qe, vid))

        # Mod 4 applies only when every non-binary linear coefficient is even.
        if !is_bin && mod(lin, 2) == 1
            discard_mod_4 = true
        end

        if mod(diag + lin, 2) == 1
            add_par!(builder_mod2, vid)
        end
    end

    con_rhs = convert(Int, con.rhs)
    mod(con_rhs, 2) == 1 && negate!(builder_mod2)

    discard_mod_4 && return builder_mod2, nothing

    builder_mod4 = XorConstraintBuilder()

    for (i, vid_i) in enumerate(con_vars)
        is_bin_i = is_binary(model.vars[vid_i])
        diag_i = convert(Int, get_quad_coeff(con.qe, vid_i, vid_i))
        lin_i = convert(Int, get_lin_coeff(con.qe, vid_i))

        if !is_bin_i
            @assert mod(lin_i, 2) == 0
            if mod(lin_i ÷ 2, 2) == 1
                add_par!(builder_mod4, vid_i)
            end
        end

        ci = is_bin_i ? diag_i + lin_i : diag_i
        if (mod(ci, 4) >> 1) == 1
            add_par!(builder_mod4, vid_i)
        end

        for j in (i + 1):lastindex(con_vars)
            vid_j = con_vars[j]
            is_bin_j = is_binary(model.vars[vid_j])
            diag_j = convert(Int, get_quad_coeff(con.qe, vid_j, vid_j))
            lin_j = convert(Int, get_lin_coeff(con.qe, vid_j))
            bilin_ij = convert(Int, get_quad_coeff(con.qe, vid_i, vid_j))

            cj = is_bin_j ? diag_j + lin_j : diag_j
            if (mod(bilin_ij, 2) == 1) ⊻ (mod(ci, 2) == 1 && mod(cj, 2) == 1)
                add_conj!(builder_mod4, vid_i, vid_j)
            end
        end
    end

    (mod(con_rhs, 4) >> 1) == 1 && negate!(builder_mod4)

    return builder_mod2, builder_mod4
end

function count_builder_occurrences!(
    var_counts::Dict{VarId, Int},
    builder::XorConstraintBuilder,
)
    for (vid, is_active) in builder.par
        is_active || continue
        var_counts[vid] = get(var_counts, vid, 0) + 1
    end

    for ((vid_i, vid_j), is_active) in builder.conj
        is_active || continue
        var_counts[vid_i] = get(var_counts, vid_i, 0) + 1
        var_counts[vid_j] = get(var_counts, vid_j, 0) + 1
    end

    return var_counts
end

function build_parity_model(model::QPModel)
    con_builders = XorConstraintBuilder[]

    for con in model.cons
        builder_mod2, builder_mod4 = get_builders(model, con)
        push!(con_builders, builder_mod2)
        builder_mod4 !== nothing && push!(con_builders, builder_mod4)
    end

    var_counts = Dict{VarId, Int}()
    for builder in con_builders
        count_builder_occurrences!(var_counts, builder)
    end

    pos_to_var_id = collect(keys(var_counts))
    sort!(pos_to_var_id; by = vid -> (-var_counts[vid], vid))
    var_id_to_pos = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var_id))

    nparity_vars = length(pos_to_var_id)
    cons = [
        build(builder, nparity_vars, var_id_to_pos) for builder in con_builders
    ]

    return ParityModel(var_id_to_pos, pos_to_var_id, cons)
end

# TODO: add postsolve stack structure here
function parity_presolve!(model::QPModel)
    propagator = PropagationManager(VarId[])

    while true
        parity_model = build_parity_model(model)
        propagator = PropagationManager(parity_model.pos_to_var_id)
        reformulate_bipartite_cons!(parity_model)
        propagate!(parity_model, propagator)

        while has_unpivoted_con(parity_model)
            if has_unpivoted_xor_con(parity_model)
                gauss_jordan_xor!(parity_model)
                propagate!(parity_model, propagator)
                substitute_parity_pivots!(parity_model)
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
