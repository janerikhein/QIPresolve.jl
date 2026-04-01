const PivotIndex = Tuple{Int, Union{Int, Nothing}}
const PivotSlot = Union{Nothing, PivotIndex}

mutable struct ParityModel
    var_id_to_pos::Dict{VarId, Int}
    pos_to_var_id::Vector{Int}
    cons::Vector{XorConstraint}
    pivots::Vector{PivotSlot}
    infeasible::Bool
end

function ParityModel(
    var_to_pos::Dict{VarId, Int},
    pos_to_var::Vector{Int},
    cons::Vector{XorConstraint},
)
    pivots = Vector{PivotSlot}(undef, length(cons))
    fill!(pivots, nothing)
    return ParityModel(var_to_pos, pos_to_var, cons, pivots, false)
end

function reformulate_bipartite_cons!(model::ParityModel)
    i = 1
    while i <= length(model.cons)
        con = ensure_updated!(model.cons[i])
        split = split_bipartite(con)

        if split === nothing
            i += 1
            continue
        end

        con1, con2 = split
        deleteat!(model.cons, i)
        insert!(model.cons, i, con2)
        insert!(model.cons, i, con1)

        deleteat!(model.pivots, i)
        insert!(model.pivots, i, nothing)
        insert!(model.pivots, i, nothing)

        i += 2
    end

    return model
end

function cleanup!(model::ParityModel)
    for i in length(model.cons):-1:1
        con = ensure_updated!(model.cons[i])
        constraint_nnz(con) == 0 || continue

        con.rhs && (model.infeasible = true)
        deleteat!(model.cons, i)
        deleteat!(model.pivots, i)
    end

    return model
end

function ensure_updated!(con::XorConstraint)
    con.meta.requires_update && update!(con)
    return con
end

constraint_nnz(con::XorConstraint) = con.meta.nnz_par + con.meta.nnz_conj

function has_unpivoted_con(model::ParityModel)
    for (i, con) in enumerate(model.cons)
        model.pivots[i] !== nothing && continue
        ensure_updated!(con)
        constraint_nnz(con) == 0 && continue
        return true
    end

    return false
end

function has_unpivoted_xor_con(model::ParityModel)
    for (i, con) in enumerate(model.cons)
        model.pivots[i] !== nothing && continue
        ensure_updated!(con)
        con.meta.is_pure_xor || continue
        constraint_nnz(con) == 0 && continue
        return true
    end

    return false
end

function is_selected_row_type(
    con::XorConstraint,
    include_xor::Bool,
    include_xor_and::Bool,
)
    ensure_updated!(con)
    return (con.meta.is_pure_xor && include_xor) || (!con.meta.is_pure_xor && include_xor_and)
end

function eliminate_pivot_from_rows!(
    model::ParityModel,
    piv_row_idx::Int,
    piv_col_idx1::Int,
    piv_col_idx2::Union{Int, Nothing},
    include_xor::Bool,
    include_xor_and::Bool,
)
    piv_con = ensure_updated!(model.cons[piv_row_idx])

    for (i, con) in enumerate(model.cons)
        i == piv_row_idx && continue
        is_selected_row_type(con, include_xor, include_xor_and) || continue
        constraint_nnz(con) == 0 && continue

        if piv_col_idx2 === nothing
            con.par[piv_col_idx1] && xor_con!(con, piv_con)
        else
            con.conj !== nothing && con.conj[piv_col_idx1, piv_col_idx2] && xor_con!(con, piv_con)
        end
    end

    return model
end

function substitute_parity_pivots!(model::ParityModel)
    for (piv_row_idx, pivot) in enumerate(model.pivots)
        pivot === nothing && continue

        piv_col_idx1, piv_col_idx2 = pivot
        piv_col_idx2 === nothing || continue

        piv_con = ensure_updated!(model.cons[piv_row_idx])
        piv_is_pure_xor = piv_con.meta.is_pure_xor

        for (row_idx, con) in enumerate(model.cons)
            row_idx == piv_row_idx && continue

            ensure_updated!(con)
            constraint_nnz(con) == 0 && continue
            con.par[piv_col_idx1] || continue

            row_pivot = model.pivots[row_idx]
            if row_pivot === nothing
                xor_con!(con, piv_con)
                continue
            end

            if con.meta.is_pure_xor == piv_is_pure_xor
                @assert !con.par[piv_col_idx1]
            else
                xor_con!(con, piv_con)
            end
        end
    end

    return model
end

"""
Applies substitution vid -> substvid ⊻ neg to all constraints
"""
function substitute_var!(model::ParityModel, vid::VarId, substid::VarId, neg::Bool)
    vid_idx = model.var_id_to_pos[vid]
    subst_idx = model.var_id_to_pos[substid]

    for con in model.cons
        con.conj === nothing && continue
        substitute_var!(con, vid_idx, subst_idx, neg)
    end

    return model
end

function get_pivot_row_idx(model::ParityModel, include_xor::Bool, include_xor_and::Bool)
    piv_row_idx = 0
    min_nnz = typemax(Int)

    for (i, con) in enumerate(model.cons)
        model.pivots[i] !== nothing && continue
        is_selected_row_type(con, include_xor, include_xor_and) || continue

        nnz = constraint_nnz(con)
        nnz == 0 && continue

        if nnz < min_nnz
            piv_row_idx = i
            min_nnz = nnz
        end
    end

    piv_row_idx == 0 && return nothing
    return piv_row_idx
end

"""
Always choose from conjunction terms if possible (i.e. not pure xor)
"""
function get_pivot_column_index(con::XorConstraint)
    ensure_updated!(con)

    if con.meta.is_pure_xor
        piv_idx = findfirst(con.par)
        @assert piv_idx !== nothing
        return piv_idx, nothing
    end

    @assert con.conj !== nothing
    conj_idx = findfirst(con.conj)
    @assert conj_idx !== nothing
    return conj_idx.I
end

function gauss_jordan!(
    model::ParityModel,
    pivot_xor::Bool,
    pivot_xor_and::Bool,
    eliminate_xor::Bool,
    eliminate_xor_and::Bool,
)
    while true
        piv_row_idx = get_pivot_row_idx(model, pivot_xor, pivot_xor_and)
        piv_row_idx === nothing && break

        piv_con = ensure_updated!(model.cons[piv_row_idx])
        piv_col_idx1, piv_col_idx2 = get_pivot_column_index(piv_con)
        model.pivots[piv_row_idx] = (piv_col_idx1, piv_col_idx2)

        eliminate_pivot_from_rows!(
            model,
            piv_row_idx,
            piv_col_idx1,
            piv_col_idx2,
            eliminate_xor,
            eliminate_xor_and,
        )
    end

    return model
end

gauss_jordan_xor!(model::ParityModel) = gauss_jordan!(model, true, false, true, false)

gauss_jordan_xor_and!(model::ParityModel) = gauss_jordan!(model, false, true, false, true)
