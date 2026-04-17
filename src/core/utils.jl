@inline _qp_iszero(val::Float64) = isapprox(val, 0.0)

@inline function _qp_num_to_str(val::Float64)
    return string(val)
end

function _qp_append_term!(io::IO, coef::Float64, body::String, is_constant::Bool, first::Bool)
    _qp_iszero(coef) && return first

    sign = coef < 0 ? "-" : "+"
    abscoef = abs(coef)

    if first
        sign == "-" && print(io, "-")
    else
        print(io, " ", sign, " ")
    end

    if is_constant
        print(io, _qp_num_to_str(abscoef))
    else
        if !isapprox(abscoef, 1.0)
            print(io, _qp_num_to_str(abscoef))
        end
        print(io, body)
    end

    return false
end

function _qp_format_expr(qe::QuadExpr)
    ids = collect(vars(qe))
    sort!(ids)

    io = IOBuffer()
    first = true
    n = length(ids)

    for i in 1:n
        id_i = ids[i]
        coef = get_quad_coeff(qe, id_i, id_i)
        if !_qp_iszero(coef)
            first = _qp_append_term!(io, coef, "x$(id_i)^2", false, first)
        end
        for j in (i + 1):n
            id_j = ids[j]
            coef = get_quad_coeff(qe, id_i, id_j)
            if !_qp_iszero(coef)
                first = _qp_append_term!(io, coef, "x$(id_i)x$(id_j)", false, first)
            end
        end
    end

    # linear terms
    for id in ids
        coef = get_lin_coeff(qe, id)
        if !_qp_iszero(coef)
            first = _qp_append_term!(io, coef, "x$(id)", false, first)
        end
    end

    # constant term
    if !_qp_iszero(qe.constant)
        first = _qp_append_term!(io, qe.constant, "", true, first)
    end

    if first
        print(io, "0")
    end

    return String(take!(io))
end

function _qp_format_constraint(con::Constraint)
    expr = _qp_format_expr(con.qe)
    lhs = con.lhs
    rhs = con.rhs

    if isfinite(lhs) && isfinite(rhs)
        return string("c", con.id, ": ", lhs, " <= ", expr, " <= ", rhs)
    elseif isfinite(lhs)
        return string("c", con.id, ": ", lhs, " <= ", expr)
    elseif isfinite(rhs)
        return string("c", con.id, ": ", expr, " <= ", rhs)
    else
        return string("c", con.id, ": ", expr)
    end
end

function _qp_show(io::IO, model::QPModel)
    println(io, "QPModel")
    println(io, "Infeasible: ", model.infeasible)

    var_ids = collect(keys(model.vars))
    sort!(var_ids)

    println(io, "Variables (", length(var_ids), "):")
    for id in var_ids
        info = model.vars[id]
        println(io, "  x", id, " in [", info.lb, ", ", info.ub, "]")
    end

    println(io, "Objective (", string(model.obj_sense), "):")
    println(io, "  ", _qp_format_expr(model.obj_expr))

    println(io, "Constraints (", length(model.cons), "):")
    for con in model.cons
        println(io, "  ", _qp_format_constraint(con))
    end
    return
end

function Base.show(io::IO, model::QPModel)
    return _qp_show(io, model)
end

function Base.show(io::IO, ::MIME"text/plain", model::QPModel)
    return _qp_show(io, model)
end

function Base.show(io::IO, qe::QuadExpr)
    return println(io, _qp_format_expr(qe))
end

function Base.show(io::IO, con::Constraint)
    return println(io, _qp_format_constraint(con))
end

function _parity_format_constraint(con::XorConstraint, pos_to_var_id::Vector{Int})
    ensure_updated!(con)

    terms = String[]
    n = length(pos_to_var_id)

    @inbounds for i in 1:n
        con.par[i] && push!(terms, "p$(pos_to_var_id[i])")
    end

    if con.conj !== nothing
        @inbounds for i in 1:(n - 1)
            for j in (i + 1):n
                con.conj[i, j] || continue
                push!(terms, "(p$(pos_to_var_id[i]) ∧ p$(pos_to_var_id[j]))")
            end
        end
    end

    lhs = isempty(terms) ? "0" : join(terms, " ⊕ ")
    return string(lhs, " = ", con.rhs ? "1" : "0")
end

function _parity_pivot_counts(model::ParityModel)
    total = 0
    xor = 0
    xor_and = 0

    for pivot in model.pivots
        pivot === nothing && continue
        total += 1

        if pivot[2] === nothing
            xor += 1
        else
            xor_and += 1
        end
    end

    return total, xor, xor_and
end

function _parity_show(io::IO, model::ParityModel)
    xor_count = 0
    xor_and_count = 0

    for con in model.cons
        ensure_updated!(con)
        if con.meta.is_pure_xor
            xor_count += 1
        else
            xor_and_count += 1
        end
    end

    pivot_total, xor_pivots, xor_and_pivots = _parity_pivot_counts(model)

    println(io, "ParityModel")
    println(io, "Variables: ", length(model.pos_to_var_id))
    println(io, "Constraints: ", length(model.cons))
    println(io, "XOR constraints: ", xor_count)
    println(io, "XOR-AND constraints: ", xor_and_count)
    println(io, "Infeasible: ", model.infeasible)
    println(io, "Pivoted rows: ", pivot_total)
    println(io, "XOR pivots: ", xor_pivots)
    println(io, "XOR-AND pivots: ", xor_and_pivots)

    println(io, "XOR constraints (", xor_count, "):")
    for con in model.cons
        ensure_updated!(con)
        con.meta.is_pure_xor || continue
        println(io, "  ", _parity_format_constraint(con, model.pos_to_var_id))
    end

    println(io, "XOR-AND constraints (", xor_and_count, "):")
    for con in model.cons
        ensure_updated!(con)
        con.meta.is_pure_xor && continue
        println(io, "  ", _parity_format_constraint(con, model.pos_to_var_id))
    end

    return
end

function Base.show(io::IO, model::ParityModel)
    return _parity_show(io, model)
end

function Base.show(io::IO, ::MIME"text/plain", model::ParityModel)
    return _parity_show(io, model)
end
