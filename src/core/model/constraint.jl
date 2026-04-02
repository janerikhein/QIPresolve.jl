"""
    Constraint

Represent a bound constraint

    lhs <= qe <= rhs

where `qe` is a `QuadExpr` and `lhs`/`rhs` are scalar bounds.
"""
mutable struct Constraint
    id::Int
    qe::QuadExpr
    lhs::Float64
    rhs::Float64
end

is_equality(con::Constraint) = (con.lhs == con.rhs)

@inline _is_integral_value(x::Float64) = isfinite(x) && x == round(x)
@inline _is_integral_or_infinite(x::Float64) = !isfinite(x) || _is_integral_value(x)
@inline _accumulate_gcd(g::Int, x::Float64) = gcd(g, abs(round(Int, x)))

function is_integer(con::Constraint)
    quad_mat = quad(con.qe)
    lin_vec = lin(con.qe)
    nvars = con.qe.nvars

    for i in 1:nvars
        _is_integral_value(lin_vec[i]) || return false

        diag_coeff = quad_mat[i, i]
        _is_integral_value(diag_coeff) || return false

        for j in (i + 1):nvars
            bilinear_coeff = quad_mat[i, j]
            _is_integral_value(bilinear_coeff) || return false
            iseven(round(Int, bilinear_coeff)) || return false
        end
    end

    lhs = con.lhs - con.qe.constant
    rhs = con.rhs - con.qe.constant
    return _is_integral_or_infinite(lhs) && _is_integral_or_infinite(rhs)
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

    quad_mat = quad(con.qe)
    lin_vec = lin(con.qe)
    nvars = con.qe.nvars
    g = 0

    for i in 1:nvars
        lin_coeff = lin_vec[i]
        lin_coeff == 0.0 || (g = _accumulate_gcd(g, lin_coeff))
        g == 1 && return false

        diag_coeff = quad_mat[i, i]
        diag_coeff == 0.0 || (g = _accumulate_gcd(g, diag_coeff))
        g == 1 && return false

        for j in (i + 1):nvars
            bilinear_coeff = quad_mat[i, j]
            bilinear_coeff == 0.0 && continue
            half_coeff = bilinear_coeff / 2.0
            _is_integral_value(half_coeff) || return false
            g = _accumulate_gcd(g, half_coeff)
            g == 1 && return false
        end
    end

    g <= 1 && return false

    scale = Float64(g)
    quad_mat ./= scale
    lin_vec ./= scale

    if isfinite(con.lhs)
        con.lhs = ceil(con.lhs / scale)
        con.lhs == 0.0 && (con.lhs = 0.0)
    end
    if isfinite(con.rhs)
        con.rhs = floor(con.rhs / scale)
        con.rhs == 0.0 && (con.rhs = 0.0)
    end

    return true
end

vars(con::Constraint) = vars(con.qe)

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
    # get views of quadratic terms and linear terms
    quad_mat = quad(con.qe)
    lin_vec = lin(con.qe)

    # symmetrize into 2*lhs <= x'(Q+Q')x + 2ax <= 2*rhs
    quad_mat .+= transpose(quad_mat)
    lin_vec .*= 2
    con.qe.constant *= 2
    con.lhs *= 2
    con.rhs *= 2
    return con
end

"""
    affine_transform!(con, var_id, scale, offset) -> Constraint

Apply the affine substitution `x_var_id := scale * x_var_id + offset`
to the constraint expression and normalize bounds.
"""
@inline function affine_transform!(con::Constraint, var_id::VarId, scale::Float64, offset::Float64)
    affine_transform!(con.qe, var_id, scale, offset)
    return normalize!(con)
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
    lin_transform!(con.qe, var_id, other_id, a, b; c = c)
    return normalize!(con)
end
