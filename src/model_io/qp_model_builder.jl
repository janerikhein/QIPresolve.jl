using ..PresolvingCore: QPModel, QuadExpr, Constraint, IntVar
import MathOptInterface as MOI

"""
    MOIConSense

Represent supported MOI constraint-set types for scalar constraints.
"""
const MOIConSense = Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo, MOI.Interval}
"""
    MOIVarSense

Represent supported MOI variable-set types for variable declarations.
"""
const MOIVarSense = Union{MOI.Integer, MOI.ZeroOne}
const VarSenseMapping = Dict(MOI.Integer => :int, MOI.ZeroOne => :bin)
const ObjSenseMapping = Dict(
    MOI.MIN_SENSE => :min,
    MOI.MAX_SENSE => :max,
    MOI.FEASIBILITY_SENSE => :feas,
)

"""
    VarId

Represent variable identifiers used by the model builders.
"""
const VarId = Int
const _LinDict = Dict{VarId, Float64}
const _QuadDict = Dict{Tuple{VarId, VarId}, Float64}

"""
    VarInfo

Store declared bounds and type information for one variable.

# Fields
- `lb`: Lower bound.
- `ub`: Upper bound.
- `var_type`: Variable type symbol, one of `:cont`, `:int`, or `:bin`.

# Throws
- `ArgumentError` if `var_type` is unsupported.
"""
struct VarInfo
    lb::Float64
    ub::Float64
    var_type::Symbol

    function VarInfo(lb::Float64, ub::Float64, var_type::Symbol)
        var_type in (:cont, :int, :bin) || throw(ArgumentError("invalid var type: $var_type"))

        return new(lb, ub, var_type)
    end
end

"""
    QuadExprBuilder

Accumulate quadratic-expression terms before materialization.

The builder stores linear and quadratic coefficients in dictionaries keyed by
variable ids and accumulates duplicate terms additively.

# Fields
- `constant`: Constant term.
- `lin`: Linear coefficients by variable id.
- `quad`: Upper-triangular quadratic coefficients by variable-id pair.
"""
mutable struct QuadExprBuilder
    constant::Float64
    lin::_LinDict
    quad::_QuadDict
end

QuadExprBuilder() = QuadExprBuilder(0, _LinDict(), _QuadDict())
QuadExprBuilder(constant::Float64) = QuadExprBuilder(constant, _LinDict(), _QuadDict())


"""
    add_lin!(expr, id, val)

Accumulate a linear coefficient into `expr`.

# Arguments
- `expr`: Quadratic-expression builder mutated in place.
- `id`: Variable id receiving the coefficient.
- `val`: Increment to add.

# Returns
- The updated coefficient stored for `id`.
"""
function add_lin!(expr::QuadExprBuilder, id::VarId, val::Float64)
    return expr.lin[id] = get(expr.lin, id, 0.0) + val
end


"""
    add_quad!(expr, id1, id2, val)

Accumulate a quadratic coefficient into `expr`.

Store the coefficient under the canonical ordered key `(min(id1, id2),
max(id1, id2))`.

# Returns
- The updated coefficient stored for the canonical variable pair.
"""
function add_quad!(expr::QuadExprBuilder, id1::VarId, id2::VarId, val::Float64)
    id1, id2 = id1 > id2 ? (id2, id1) : (id1, id2)
    return expr.quad[(id1, id2)] = get(expr.quad, (id1, id2), 0) + val
end


"""
    build(builder::QuadExprBuilder)

Materialize `builder` as a `QuadExpr`.

Drop near-zero stored coefficients and construct a `PresolvingCore.QuadExpr`
from the remaining terms.

# Returns
- A `QuadExpr` equivalent to the accumulated builder contents.
"""
function build(builder::QuadExprBuilder)
    quad_terms = Vector{Tuple{Float64, VarId, VarId}}()
    sizehint!(quad_terms, length(builder.quad))
    @inbounds for ((id1, id2), val) in builder.quad
        isapprox(val, 0.0) && continue
        push!(quad_terms, (val, id1, id2))
    end

    lin_terms = Vector{Tuple{Float64, VarId}}()
    sizehint!(lin_terms, length(builder.lin))
    @inbounds for (id, val) in builder.lin
        isapprox(val, 0.0) && continue
        push!(lin_terms, (val, id))
    end

    return QuadExpr(quad_terms, lin_terms; constant = builder.constant)
end


"""
    ConstraintBuilder

Store one not-yet-materialized quadratic constraint.

# Fields
- `id`: Constraint identifier.
- `expr`: Expression builder for the left-hand side.
- `lhs`: Lower bound.
- `rhs`: Upper bound.
"""
struct ConstraintBuilder
    id::Int
    expr::QuadExprBuilder
    lhs::Float64
    rhs::Float64
end

"""
    build(builder::ConstraintBuilder)

Materialize `builder` as a `Constraint`.

# Returns
- A `PresolvingCore.Constraint` built from the stored expression and bounds.
"""
function build(builder::ConstraintBuilder)
    quad_expr = build(builder.expr)

    return Constraint(builder.id, quad_expr, builder.lhs, builder.rhs)
end

"""
    QPModelBuilder

Accumulate a quadratic-program model before final construction.

# Fields
- `vars`: Declared variable metadata keyed by variable id.
- `cons`: Constraint builders in insertion order.
- `obj_expr`: Objective-expression builder.
- `obj_sense`: Objective sense symbol.
- `infeasible`: Whether imported metadata already proves infeasibility.
"""
mutable struct QPModelBuilder
    vars::Dict{VarId, VarInfo}
    cons::Vector{ConstraintBuilder}
    obj_expr::QuadExprBuilder
    obj_sense::Symbol
    infeasible::Bool
end


QPModelBuilder() = QPModelBuilder(
    Dict{VarId, VarInfo}(),
    Vector{ConstraintBuilder}(),
    QuadExprBuilder(),
    :undef,
    false,
)

"""
    from_moi(moi_model)

Translate an MOI model into a `QPModelBuilder`.

Import the objective and every supported constraint family from `moi_model`
into the builder representation used by QIPresolve.

# Returns
- A populated `QPModelBuilder`.
"""
function from_moi(moi_model)
    qp_model = QPModelBuilder()

    # register objective
    obj_expr_type = MOI.get(moi_model, MOI.ObjectiveFunctionType())
    obj_expr = MOI.get(moi_model, MOI.ObjectiveFunction{obj_expr_type}())
    obj_sense::MOI.OptimizationSense = MOI.get(moi_model, MOI.ObjectiveSense())
    register_objective!(qp_model, obj_expr, obj_sense)

    # register constraints
    for (expr_type, sense_type) in MOI.get(moi_model, MOI.ListOfConstraintTypesPresent())
        con_indices = MOI.get(moi_model, MOI.ListOfConstraintIndices{expr_type, sense_type}())
        for ci in con_indices
            f = MOI.get(moi_model, MOI.ConstraintFunction(), ci)::expr_type
            s = MOI.get(moi_model, MOI.ConstraintSet(), ci)::sense_type
            register_constraint!(qp_model, f, s)
        end
    end

    return qp_model
end

"""
    register_con!(model; quad_expr_terms=[], lin_expr_terms=[], constant=0.0, lhs=-Inf, rhs=Inf)

Register a new constraint in `model`.

Accumulate the supplied linear and quadratic terms, ensure their variables are
present in `model.vars`, and append a new `ConstraintBuilder`.

# Returns
- The newly appended `ConstraintBuilder`.

# Side Effects
- Mutates `model.vars` and `model.cons`.
"""
function register_con!(
        model::QPModelBuilder;
        quad_expr_terms::Vector{Tuple{Float64, VarId, VarId}} = Tuple{Float64, VarId, VarId}[],
        lin_expr_terms::Vector{Tuple{Float64, VarId}} = Tuple{Float64, VarId}[],
        constant::Float64 = 0.0,
        lhs::Float64 = -Inf,
        rhs::Float64 = Inf
    )
    quad_expr = QuadExprBuilder(constant)

    for (coeff, var_id_1, var_id_2) in quad_expr_terms
        register_var_info!(model, var_id_1)
        register_var_info!(model, var_id_2)
        add_quad!(quad_expr, var_id_1, var_id_2, coeff)
    end

    for (coeff, var_id) in lin_expr_terms
        register_var_info!(model, var_id)
        add_lin!(quad_expr, var_id, coeff)
    end

    con_id = length(model.cons) + 1

    return push!(model.cons, ConstraintBuilder(con_id, quad_expr, lhs, rhs))
end

"""
    register_var_info!(model, id; var_type=:cont, lb=-Inf, ub=Inf)

Register or tighten metadata for variable `id`.

When `id` already exists, the stored bounds are intersected with `[lb, ub]` and
the variable type is restricted to the stronger existing type.

# Returns
- The updated `VarInfo` stored for `id`.
"""
function register_var_info!(
        model::QPModelBuilder, id::VarId;
        var_type::Symbol = :cont, lb::Float64 = -Inf, ub::Float64 = Inf
    )
    if !haskey(model.vars, id)
        lb > ub && (model.infeasible = true)
        model.vars[id] = VarInfo(lb, ub, var_type)
        return
    else
        var_info = model.vars[id]
    end

    # restrict var_type
    if (var_info.var_type == :bin && (var_type == :int || var_type == :cont)) || (var_info.var_type == :int && var_type == :cont)
        var_type = var_info.var_type
    end

    # restrict bounds
    lb = max(var_info.lb, lb)
    ub = min(var_info.ub, ub)
    lb > ub && (model.infeasible = true)

    return model.vars[id] = VarInfo(lb, ub, var_type)
end

"""
    register_obj!(model, constant, obj_sense; quad_expr_terms=[], lin_expr_terms=[])

Register the objective function for `model`.

Replace the current objective builder with one built from the supplied terms and
store the requested objective sense.

# Returns
- The updated objective-sense symbol.

# Side Effects
- Mutates `model.obj_expr`, `model.obj_sense`, and may register new variables.
"""
function register_obj!(
        model::QPModelBuilder,
        constant::Float64,
        obj_sense::Symbol;
        quad_expr_terms::Vector{Tuple{Float64, VarId, VarId}} = Tuple{Float64, VarId, VarId}[],
        lin_expr_terms::Vector{Tuple{Float64, VarId}} = Tuple{Float64, VarId}[],
    )
    quad_expr = QuadExprBuilder(constant)

    for (coeff, var_id_1, var_id_2) in quad_expr_terms
        register_var_info!(model, var_id_1)
        register_var_info!(model, var_id_2)
        add_quad!(quad_expr, var_id_1, var_id_2, coeff)
    end

    for (coeff, var_id) in lin_expr_terms
        register_var_info!(model, var_id)
        add_lin!(quad_expr, var_id, coeff)
    end

    model.obj_expr = quad_expr
    return model.obj_sense = obj_sense
end

"""
    get_term(moi_term)

Convert one MOI scalar term into the builder tuple format.

The affine overload returns `(coefficient, variable_id)`. The quadratic
overload returns `(coefficient, variable_id_1, variable_id_2)` and rescales
diagonal terms to undo MOF's doubled-diagonal convention.
"""
get_bound(moi_sense::MOI.LessThan{T}) where {T <: Real} = (-Inf, Float64(moi_sense.upper))
get_bound(moi_sense::MOI.GreaterThan{T}) where {T <: Real} = (Float64(moi_sense.lower), Inf)
get_bound(moi_sense::MOI.EqualTo{T}) where {T <: Real} = (moi_sense.value, moi_sense.value) .|> Float64
get_bound(moi_sense::MOI.Interval{T}) where {T <: Real} = (moi_sense.lower, moi_sense.upper) .|> Float64

get_term(moi_term::MOI.ScalarAffineTerm{T}) where {T <: Real} = (Float64(moi_term.coefficient), Int(moi_term.variable.value))

# note: somehow mof files scale quad diagonal entries by 2
function get_term(moi_term::MOI.ScalarQuadraticTerm{T}) where {T <: Real}
    coeff = Float64(moi_term.coefficient)
    var1 = Int(moi_term.variable_1.value)
    var2 = Int(moi_term.variable_2.value)
    if var1 == var2
        coeff = coeff / 2
    end

    return (coeff, var1, var2)
end

"""
    _parse_quad_expr(expr)

Extract quadratic and linear term tuples from an MOI quadratic expression.

# Returns
- A tuple `(quad_terms, lin_terms)` in `QPModelBuilder` input format.
"""
function _parse_quad_expr(expr::MOI.ScalarQuadraticFunction{T}) where {T <: Real}
    quad_terms = Vector{Tuple{Float64, Int, Int}}(undef, length(expr.quadratic_terms))
    @inbounds for (i, term) in enumerate(expr.quadratic_terms)
        quad_terms[i] = get_term(term)
    end

    lin_terms = Vector{Tuple{Float64, Int}}(undef, length(expr.affine_terms))
    @inbounds for (i, term) in enumerate(expr.affine_terms)
        lin_terms[i] = get_term(term)
    end

    return quad_terms, lin_terms
end

"""
    _parse_affine_expr(expr)

Extract linear term tuples from an MOI affine expression.
"""
function _parse_affine_expr(expr::MOI.ScalarAffineFunction{T}) where {T <: Real}
    lin_terms = Vector{Tuple{Float64, Int}}(undef, length(expr.terms))
    @inbounds for (i, term) in enumerate(expr.terms)
        lin_terms[i] = get_term(term)
    end

    return lin_terms
end

"""
    register_constraint!(qp_model, expr, sense)

Register one MOI constraint in `qp_model`.

Supported overloads handle affine scalar functions, quadratic scalar
functions, scalar variable bounds, and integer/binary variable declarations.
"""
function register_constraint!(qp_model::QPModelBuilder, expr::MOI.ScalarAffineFunction{T}, sense::MOIConSense) where {T <: Real}
    lin_terms = _parse_affine_expr(expr)
    lhs, rhs = get_bound(sense)
    constant = expr.constant

    return register_con!(qp_model; lin_expr_terms = lin_terms, constant = constant, lhs = lhs, rhs = rhs)
end


function register_constraint!(qp_model::QPModelBuilder, expr::MOI.ScalarQuadraticFunction{T}, sense::MOIConSense) where {T <: Real}
    quad_terms, lin_terms = _parse_quad_expr(expr)
    lhs, rhs = get_bound(sense)
    constant = expr.constant

    return register_con!(qp_model; lin_expr_terms = lin_terms, quad_expr_terms = quad_terms, constant = constant, lhs = lhs, rhs = rhs)
end


function register_constraint!(qp_model::QPModelBuilder, expr::MOI.VariableIndex, sense::MOIConSense)
    lhs, rhs = get_bound(sense)

    return register_var_info!(qp_model, expr.value; lb = lhs, ub = rhs)
end


register_constraint!(qp_model::QPModelBuilder, expr::MOI.VariableIndex, ::S) where {S <: MOIVarSense} =
    register_var_info!(qp_model, expr.value; var_type = VarSenseMapping[S])


"""
    register_objective!(qp_model, expr, sense)

Register the MOI objective in `qp_model`.

Convert the supported affine objective into builder form and map the MOI
optimization sense to the corresponding QIPresolve symbol.
"""
function register_objective!(qp_model::QPModelBuilder, expr::MOI.ScalarAffineFunction{T}, sense::MOI.OptimizationSense) where {T <: Real}
    lin_terms = _parse_affine_expr(expr)
    constant = expr.constant
    obj_sense = ObjSenseMapping[sense]

    return register_obj!(qp_model, constant, obj_sense; lin_expr_terms = lin_terms)
end


function register_objective!(qp_model::QPModelBuilder, expr::MOI.ScalarQuadraticFunction{T}, sense::MOI.OptimizationSense) where {T <: Real}
    quad_terms, lin_terms = _parse_quad_expr(expr)
    constant = expr.constant
    obj_sense = ObjSenseMapping[sense]

    return register_obj!(qp_model, constant, obj_sense; quad_expr_terms = quad_terms, lin_expr_terms = lin_terms)
end


"""
    build_model(builder)

Materialize `builder` as a `QPModel`.

Convert builder-side variable metadata and constraints into
`PresolvingCore.IntVar`, `Constraint`, and `QuadExpr` objects.

# Returns
- A `PresolvingCore.QPModel`.

# Throws
- `ErrorException` if a continuous variable is present, because the presolve
  core currently supports only integer and binary variables.
"""
function build_model(builder::QPModelBuilder)

    vars = Dict{VarId, IntVar}()
    sizehint!(vars, length(builder.vars))
    for (var_id, var_info) in builder.vars
        var_info.var_type == :cont && error("continous variables are unsupported")
        vars[var_id] = IntVar(var_info.lb, var_info.ub)
    end

    cons = Vector{Constraint}(undef, length(builder.cons))
    @inbounds for (i, con_builder) in enumerate(builder.cons)
        cons[i] = build(con_builder)
    end

    obj_expr = build(builder.obj_expr)
    obj_sense = builder.obj_sense

    model = QPModel(vars, cons, obj_expr, obj_sense)
    model.infeasible = builder.infeasible
    return model
end
