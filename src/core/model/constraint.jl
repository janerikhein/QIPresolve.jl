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

    function Constraint(id::Int, qe::QuadExpr, lhs::Float64, rhs::Float64)
        con = new(id, qe, lhs, rhs)
        normalize!(con)
        return con
    end
end

"""
    is_equality(con)

Check whether `con` has identical lower and upper bounds.
"""
is_equality(con::Constraint) = (con.lhs == con.rhs)

"""
    is_integer(con)

Check whether every coefficient in `con.qe` satisfies the integral-parity
invariants required by the presolve routines.
"""
@inline is_integer(con::Constraint) = isinteger(con.qe)

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
        if lin_coeff != 0.0
            g = gcd(g, abs(trunc(Int, lin_coeff)))
        end
        g == 1 && return false

        diag_coeff = get_quad_coeff(con.qe, vid_i, vid_i)
        if diag_coeff != 0.0
            g = gcd(g, abs(trunc(Int, diag_coeff)))
        end
        g == 1 && return false

        for j in (i + 1):nvars
            vid_j = var_ids[j]
            bilinear_coeff = get_quad_coeff(con.qe, vid_i, vid_j)
            bilinear_coeff == 0.0 && continue
            g = gcd(g, abs(trunc(Int, bilinear_coeff) ÷ 2))
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
        con.lhs = _canonicalize_zero(ceil(con.lhs / scale))
    end
    if isfinite(con.rhs)
        con.rhs = _canonicalize_zero(floor(con.rhs / scale))
    end

    return true
end

vars(con::Constraint) = vars(con.qe)

"""
    remove_var!(con, id; clear_buf=true)

Remove variable `id` from the quadratic expression of `con`.

# Returns
- `true` if `id` was present and removed.
- `false` if `id` was absent.

# Side Effects
- Mutates `con.qe`.
"""
function remove_var!(con::Constraint, id::VarId; clear_buf::Bool = true)
    removed = remove_var!(con.qe, id; clear_buf = clear_buf)
    removed || return false
    return true
end

"""
    normalize!(con::Constraint) -> Constraint

Move the constant term from the quadratic expression into the bounds.

After normalization, `con.qe.constant == 0.0` and `lhs`/`rhs` are shifted by the
previous constant so the constraint is equivalent. Integral expressions also
tighten finite bounds to integers.
"""
function normalize!(con::Constraint)
    normalize!(con.qe)
    con.lhs -= con.qe.constant
    con.rhs -= con.qe.constant
    con.qe.constant = 0.0

    if isinteger(con.qe)
        isfinite(con.lhs) && (con.lhs = ceil(con.lhs))
        isfinite(con.rhs) && (con.rhs = floor(con.rhs))
    end

    con.lhs = _canonicalize_zero(con.lhs)
    con.rhs = _canonicalize_zero(con.rhs)
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
    quad_mat = quad(con.qe)
    lin_vec = lin(con.qe)

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
    has_var(con.qe, var_id) || return con

    if isinteger(con.qe)
        isinteger(scale) || throw(ArgumentError("Integral constraints require integer-valued scale, got $scale"))
        isinteger(offset) || throw(ArgumentError("Integral constraints require integer-valued offset, got $offset"))
    end

    affine_transform!(con.qe, var_id, scale, offset)
    normalize!(con)
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

    if isinteger(con.qe)
        isinteger(a) || throw(ArgumentError("Integral constraints require integer-valued a, got $a"))
        isinteger(b) || throw(ArgumentError("Integral constraints require integer-valued b, got $b"))
        isinteger(c) || throw(ArgumentError("Integral constraints require integer-valued c, got $c"))
    end

    if b != 0.0 && var_id != other_id && has_var(con.qe, var_id) && !has_var(con.qe, other_id)
        add_var!(con.qe, other_id; clear_buf = true)
    end
    lin_transform!(con.qe, var_id, other_id, a, b; c = c)
    normalize!(con)
    return con
end
