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
    lin_transform!(con, var_id, other_id, a, b) -> Constraint

Apply the linear substitution `x_var_id := a * x_var_id + b * x_other_id`
to the constraint expression and normalize bounds.
"""
@inline function lin_transform!(con::Constraint, var_id::VarId, other_id::VarId, a::Float64, b::Float64)
    lin_transform!(con.qe, var_id, other_id, a, b)
    return normalize!(con)
end
