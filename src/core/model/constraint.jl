"""
    Constraint

Represent a bound constraint

    lhs <= qe <= rhs

where `qe` is a `QuadExpr` and `lhs`/`rhs` are scalar bounds.
"""
@inline _is_integral_value(x::Float64) = isfinite(x) && x == round(x)
@inline _is_integral_or_infinite(x::Float64) = !isfinite(x) || _is_integral_value(x)
@inline _accumulate_gcd(g::Int, x::Float64) = gcd(g, abs(round(Int, x)))

function _compute_is_integer(qe::QuadExpr, lhs::Float64, rhs::Float64)
    var_ids = vars(qe)
    nvars = length(var_ids)

    for i in 1:nvars
        vid_i = var_ids[i]

        lin_coeff = get_lin_coeff(qe, vid_i)
        _is_integral_value(lin_coeff) || return false

        diag_coeff = get_quad_coeff(qe, vid_i, vid_i)
        _is_integral_value(diag_coeff) || return false

        for j in (i + 1):nvars
            vid_j = var_ids[j]
            bilinear_coeff = get_quad_coeff(qe, vid_i, vid_j)
            _is_integral_value(bilinear_coeff) || return false
            iseven(round(Int, bilinear_coeff)) || return false
        end
    end

    normalized_lhs = lhs - qe.constant
    normalized_rhs = rhs - qe.constant
    return _is_integral_or_infinite(normalized_lhs) && _is_integral_or_infinite(normalized_rhs)
end

mutable struct Constraint
    id::Int
    qe::QuadExpr
    lhs::Float64
    rhs::Float64
    is_integer::Bool

    function Constraint(id::Int, qe::QuadExpr, lhs::Float64, rhs::Float64)
        return new(id, qe, lhs, rhs, _compute_is_integer(qe, lhs, rhs))
    end
end

is_equality(con::Constraint) = (con.lhs == con.rhs)

@inline is_integer(con::Constraint) = con.is_integer

@inline function _refresh_is_integer!(con::Constraint)
    con.is_integer = _compute_is_integer(con.qe, con.lhs, con.rhs)
    return con
end

@inline function _require_integral_transform_value(name::AbstractString, value::Float64)
    _is_integral_value(value) && return
    throw(ArgumentError("Integral constraints require integer-valued $name, got $value"))
end

@inline function _require_integer_preserving_affine(scale::Float64, offset::Float64)
    _require_integral_transform_value("scale", scale)
    _require_integral_transform_value("offset", offset)
    return
end

@inline function _require_integer_preserving_lin(a::Float64, b::Float64, c::Float64)
    _require_integral_transform_value("a", a)
    _require_integral_transform_value("b", b)
    _require_integral_transform_value("c", c)
    return
end

"""
    scale_gcd!(con::Constraint) -> Bool

Normalize `con`, divide all coefficients by their greatest common divisor, and
tighten finite bounds to integers. Returns `true` iff the constraint was scaled.
Constraints with non-integral data are left unchanged after normalization.
"""
function scale_gcd!(con::Constraint)
    normalize!(con)
    is_integer(con) || return false

    var_ids = collect(vars(con.qe))
    nvars = length(var_ids)
    g = 0

    for i in 1:nvars
        vid_i = var_ids[i]

        lin_coeff = get_lin_coeff(con.qe, vid_i)
        lin_coeff == 0.0 || (g = _accumulate_gcd(g, lin_coeff))
        g == 1 && return false

        diag_coeff = get_quad_coeff(con.qe, vid_i, vid_i)
        diag_coeff == 0.0 || (g = _accumulate_gcd(g, diag_coeff))
        g == 1 && return false

        for j in (i + 1):nvars
            vid_j = var_ids[j]
            bilinear_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            bilinear_coeff == 0.0 && continue
            half_coeff = bilinear_coeff / 2.0
            _is_integral_value(half_coeff) || return false
            g = _accumulate_gcd(g, half_coeff)
            g == 1 && return false
        end
    end

    g <= 1 && return false

    scale = Float64(g)

    for vid in var_ids
        lin_coeff = get_lin_coeff(con.qe, vid)
        lin_coeff == 0.0 || set_lin_coeff!(con.qe, vid, lin_coeff / scale)

        diag_coeff = get_quad_coeff(con.qe, vid, vid)
        diag_coeff == 0.0 || set_quad_coeff!(con.qe, vid, vid, diag_coeff / scale)
    end

    for i in 1:nvars
        vid_i = var_ids[i]
        for j in (i + 1):nvars
            vid_j = var_ids[j]
            bilinear_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            bilinear_coeff == 0.0 || set_quad_coeff!(con.qe, vid_i, vid_j, bilinear_coeff / scale)
        end
    end

    if isfinite(con.lhs)
        con.lhs = ceil(con.lhs / scale)
        con.lhs == 0.0 && (con.lhs = 0.0)
    end
    if isfinite(con.rhs)
        con.rhs = floor(con.rhs / scale)
        con.rhs == 0.0 && (con.rhs = 0.0)
    end

    _refresh_is_integer!(con)
    return true
end

vars(con::Constraint) = vars(con.qe)

function remove_var!(con::Constraint, id::VarId; clear_buf::Bool = true)
    removed = remove_var!(con.qe, id; clear_buf = clear_buf)
    removed || return false
    con.is_integer || _refresh_is_integer!(con)
    return true
end

"""
    normalize!(con::Constraint) -> Constraint

Move the constant term from the quadratic expression into the bounds.

After normalization, `con.qe.constant == 0.0` and `lhs`/`rhs` are shifted by the
previous constant so the constraint is equivalent.
"""
function normalize!(con::Constraint)
    normalize!(con.qe)
    con.lhs -= con.qe.constant
    con.rhs -= con.qe.constant
    con.qe.constant = 0.0
    return con
end


"""
    symmetrize!(con::Constraint) -> Constraint

Symmetrize the quadratic form in-place.

For the dense-buffer representation, this replaces `Q` by `Q + Q'` and rescales
the linear term and bounds so the constraint remains equivalent:

    lhs <= x'Qx + c'x + constant <= rhs

becomes

    2lhs <= x'(Q+Q')x + 2c'x + 2constant <= 2rhs
"""
function symmetrize!(con::Constraint)
    was_integer = con.is_integer

    # get views of quadratic terms and linear terms
    quad_mat = quad(con.qe)
    lin_vec = lin(con.qe)

    # symmetrize into 2*lhs <= x'(Q+Q')x + 2ax <= 2*rhs
    quad_mat .+= transpose(quad_mat)
    lin_vec .*= 2
    con.qe.constant *= 2
    con.lhs *= 2
    con.rhs *= 2
    was_integer || _refresh_is_integer!(con)
    return con
end

"""
    affine_transform!(con, var_id, scale, offset) -> Constraint

Apply the affine substitution `x_var_id := scale * x_var_id + offset`
to the constraint expression and normalize bounds.
"""
@inline function affine_transform!(con::Constraint, var_id::VarId, scale::Float64, offset::Float64)
    has_var(con.qe, var_id) || return con

    was_integer = con.is_integer
    was_integer && _require_integer_preserving_affine(scale, offset)

    affine_transform!(con.qe, var_id, scale, offset)
    normalize!(con)
    was_integer || _refresh_is_integer!(con)
    return con
end

"""
    lin_transform!(con, var_id, other_id, a, b; c = 0.0) -> Constraint

Apply the linear substitution `x_var_id := a * x_var_id + b * x_other_id + c`
to the constraint expression and normalize bounds.
"""
@inline function lin_transform!(
    con::Constraint,
    var_id::VarId,
    other_id::VarId,
    a::Float64,
    b::Float64;
    c::Float64 = 0.0,
)
    has_var(con.qe, var_id) || return con
    var_id == other_id && return con

    was_integer = con.is_integer
    was_integer && _require_integer_preserving_lin(a, b, c)

    if b != 0.0 && var_id != other_id && has_var(con.qe, var_id) && !has_var(con.qe, other_id)
        add_var!(con.qe, other_id; clear_buf = true)
    end
    lin_transform!(con.qe, var_id, other_id, a, b; c = c)
    normalize!(con)
    was_integer || _refresh_is_integer!(con)
    return con
end
