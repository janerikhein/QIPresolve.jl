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

function substitute_pivots_in_conjunctive_terms!(model::ParityModel)
    changed = falses(length(model.cons))

    for (piv_row_idx, pivot) in enumerate(model.pivots)
        pivot === nothing && continue

        piv_col_idx1, piv_col_idx2 = pivot
        piv_col_idx2 === nothing || continue

        piv_con = ensure_updated!(model.cons[piv_row_idx])
        @assert piv_con.meta.is_pure_xor
        @assert piv_con.par[piv_col_idx1]

        subst_mask = copy(piv_con.par)
        subst_mask[piv_col_idx1] = false

        for (row_idx, con) in enumerate(model.cons)
            row_idx == piv_row_idx && continue

            ensure_updated!(con)
            constraint_nnz(con) == 0 && continue
            con.conj === nothing && continue
            any(@view con.conj[:, piv_col_idx1]) || continue

            substitute_var_in_conjunctive_terms!(con, piv_col_idx1, subst_mask, piv_con.rhs)
            changed[row_idx] = true
        end
    end

    for row_idx in eachindex(model.cons)
        changed[row_idx] || continue
        _revalidate_pivot!(model, row_idx)
    end

    return model
end

function _constraint_contains_var(con::XorConstraint, var_idx::Int)
    con.par[var_idx] && return true
    con.conj === nothing && return false
    return any(@view con.conj[:, var_idx])
end

function _fix_var_rows!(changed::BitVector, model::ParityModel, vid::VarId, val::Bool)
    
    vid_idx = model.var_id_to_pos[vid]

    for (i, con) in enumerate(model.cons)
        _constraint_contains_var(con, vid_idx) || continue
        fix_var!(con, vid_idx, val)
        changed[i] = true
    end

    return model
end

function _substitute_var_rows!(
    changed::BitVector,
    model::ParityModel,
    vid::VarId,
    substid::VarId,
    neg::Bool,
)
    if vid == substid
        @assert !neg
        return model
    end

    vid_idx = model.var_id_to_pos[vid]
    subst_idx = model.var_id_to_pos[substid]

    for (i, con) in enumerate(model.cons)
        _constraint_contains_var(con, vid_idx) || continue
        substitute_var!(con, vid_idx, subst_idx, neg)
        changed[i] = true
    end

    return model
end

function _is_valid_pivot(con::XorConstraint, pivot::PivotIndex)
    piv_col_idx1, piv_col_idx2 = pivot
    if piv_col_idx2 === nothing
        return con.par[piv_col_idx1]
    end

    return con.conj !== nothing && con.conj[piv_col_idx1, piv_col_idx2]
end

function _revalidate_pivot!(model::ParityModel, row_idx::Int)
    pivot = model.pivots[row_idx]
    pivot === nothing && return model

    con = ensure_updated!(model.cons[row_idx])
    _is_valid_pivot(con, pivot) || (model.pivots[row_idx] = nothing)

    return model
end

function _has_constraints_requiring_propagation(model::ParityModel)
    return any(con.meta.requires_prop for con in model.cons)
end

_assignment_lit(vid::VarId, val::Bool) = VarLit(vid, !val)

function _add_tripartite_implications!(
    manager::PropagationManager,
    idx_to_vid::Vector{VarId},
    action::_TripartiteImplicationAction,
)
    x_vid = idx_to_vid[action.x_idx]
    y_vid = idx_to_vid[action.y_idx]

    add_implication!(
        manager,
        _assignment_lit(x_vid, !action.r1),
        _assignment_lit(y_vid, action.r2),
    )
    add_implication!(
        manager,
        _assignment_lit(y_vid, !action.r2),
        _assignment_lit(x_vid, action.r1),
    )

    return manager
end

function _insert_tripartite_rewrite!(
    model::ParityModel,
    manager::PropagationManager,
    row_idx::Int,
    action::_TripartiteRewriteAction,
)
    deleteat!(model.cons, row_idx)
    deleteat!(model.pivots, row_idx)

    new_cons = XorConstraint[]
    for (support, rhs) in (
        (action.support1, action.rhs1),
        (action.support2, action.rhs2),
    )
        nnz = count(support)
        @assert nnz > 0

        if nnz == 1
            vid = model.pos_to_var_id[findfirst(support)]
            fix_var!(manager, vid, rhs)
        else
            push!(new_cons, XorConstraint(copy(support), rhs))
        end
    end

    for (offset, con) in enumerate(new_cons)
        insert!(model.cons, row_idx + offset - 1, con)
        insert!(model.pivots, row_idx + offset - 1, nothing)
    end

    return model
end

"""
Applies fixing vid -> val to all constraints
"""
function fix_var!(model::ParityModel, vid::VarId, val::Bool)
    changed = falses(length(model.cons))
    return _fix_var_rows!(changed, model, vid, val)
end

"""
Applies substitution vid -> substvid ⊻ neg to all constraints
"""
function substitute_var!(model::ParityModel, vid::VarId, substid::VarId, neg::Bool)
    changed = falses(length(model.cons))
    return _substitute_var_rows!(changed, model, vid, substid, neg)
end

function propagate!(model::ParityModel, manager::PropagationManager)
    while _has_constraints_requiring_propagation(model) && !model.infeasible
        empty!(manager.seen_fixings)
        empty!(manager.seen_substitutions)

        i = 1
        while i <= length(model.cons)
            con = ensure_updated!(model.cons[i])
            if !con.meta.requires_prop
                i += 1
                continue
            end

            if con.meta.is_pure_xor || (con.meta.nnz_par == 0 && con.meta.nnz_conj == 1)
                register_implications!(manager, con, model.pos_to_var_id)
                i += 1
                continue
            end

            tripartite_action = _tripartite_action(con)
            if tripartite_action isa _TripartiteImplicationAction
                _add_tripartite_implications!(manager, model.pos_to_var_id, tripartite_action)
                con.meta.requires_prop = false
                i += 1
                continue
            elseif tripartite_action isa _TripartiteRewriteAction
                _insert_tripartite_rewrite!(model, manager, i, tripartite_action)
                continue
            end

            register_implications!(manager, con, model.pos_to_var_id)
            i += 1
        end

        update!(manager)
        changed = falses(length(model.cons))

        while true
            fixing = pop_fixing!(manager)
            fixing === nothing && break
            vid, val = fixing
            _fix_var_rows!(changed, model, vid, val)
        end
        
        while true
            substitution = pop_substitution!(manager)
            substitution === nothing && break
            vid, substid, neg = substitution
            _substitute_var_rows!(changed, model, vid, substid, neg)
        end

        for i in eachindex(model.cons)
            changed[i] || continue
            _revalidate_pivot!(model, i)
        end

        cleanup!(model)
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
