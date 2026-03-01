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


function Base.show(io::IO, con::XorConstraint)
    n = length(con.pos_to_var)
    terms = String[]

    # Linear XOR terms
    @inbounds for i in 1:n
        if con.par[i]
            push!(terms, "p$(con.pos_to_var[i])")
        end
    end
    if con.conj !== nothing
        # Conjunction (quadratic) terms
        @inbounds for i in 1:(n - 1)
            for j in (i + 1):n
                if con.conj[i, j]
                    v1 = con.pos_to_var[i]
                    v2 = con.pos_to_var[j]
                    push!(terms, "(p$(v1) ∧ p$(v2))")
                end
            end
        end
    end

    if isempty(terms)
        print(io, con.rhs ? "0 = 1" : "0 = 0")
        return
    end

    print(io, join(terms, " ⊕  "))
    return print(io, " = ", con.rhs ? "1" : "0")
end
