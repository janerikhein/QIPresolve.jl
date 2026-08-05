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
