using Test
using QIPresolve
import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
import QIPresolve.PresolvingCore as PC

const IOMOI = QIPresolve.ModelIO.MOI
const IOQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const IOLinTerm = Tuple{Float64, PC.VarId}

io_empty_qe() = PC.QuadExpr(IOQuadTerm[], IOLinTerm[])

function io_load_lp_core_model(file_path::AbstractString)
    moi_model = FF.Model(format = FF.FORMAT_LP)
    MOI.read_from_file(moi_model, file_path)
    return QIPresolve.build_model(QIPresolve.from_moi(moi_model))
end

function io_constraint_indices(model, ::Type{F}, ::Type{S}) where {F, S}
    return IOMOI.get(model, IOMOI.ListOfConstraintIndices{F, S}())
end

function io_assert_same_expr(actual::PC.QuadExpr, expected::PC.QuadExpr)
    actual_vars = sort!(collect(PC.vars(actual)))
    expected_vars = sort!(collect(PC.vars(expected)))
    @test actual_vars == expected_vars
    @test isapprox(actual.constant, expected.constant; atol = 1.0e-12)

    for var_id in expected_vars
        @test isapprox(PC.get_lin_coeff(actual, var_id), PC.get_lin_coeff(expected, var_id); atol = 1.0e-12)
    end
    for (i, var_id_1) in enumerate(expected_vars)
        for var_id_2 in @view expected_vars[i:end]
            @test isapprox(PC.get_quad_coeff(actual, var_id_1, var_id_2), PC.get_quad_coeff(expected, var_id_1, var_id_2); atol = 1.0e-12)
        end
    end
end

@testset "LP file import path supports presolve" begin
    lp_path = tempname() * ".lp"
    write(
        lp_path,
        """
        Minimize
         obj: x1 + 2 x2
        Subject To
         c1: x1 + x2 <= 1
        Bounds
         0 <= x1 <= 1
         0 <= x2 <= 1
        Binary
         x1
         x2
        End
        """,
    )

    try
        model = io_load_lp_core_model(lp_path)
        before_nvars = length(model.vars)
        before_ncons = length(model.cons)

        result = QIPresolve.presolve!(model)

        @test result isa QIPresolve.PresolveResult
        @test result.model === model
        @test !model.infeasible
        @test before_nvars == 2
        @test before_ncons == 1
    finally
        rm(lp_path; force = true)
    end
end

@testset "build_moi_model exports variable domains, bounds, and names" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        3 => PC.IntVar(0.0, 1.0),
        7 => PC.IntVar(-Inf, 5.0),
        9 => PC.IntVar(-2.0, Inf),
        11 => PC.IntVar(4.0, 4.0),
    )
    model = PC.QPModel(vars, PC.Constraint[], io_empty_qe(), :feas)

    moi_model = QIPresolve.build_moi_model(model)
    moi_vars = IOMOI.get(moi_model, IOMOI.ListOfVariableIndices())

    @test moi_vars == IOMOI.VariableIndex.(1:4)
    @test [IOMOI.get(moi_model, IOMOI.VariableName(), v) for v in moi_vars] == ["x3", "x7", "x9", "x11"]

    zero_one = io_constraint_indices(moi_model, IOMOI.VariableIndex, IOMOI.ZeroOne)
    integers = io_constraint_indices(moi_model, IOMOI.VariableIndex, IOMOI.Integer)
    intervals = io_constraint_indices(moi_model, IOMOI.VariableIndex, IOMOI.Interval{Float64})
    lower_bounds = io_constraint_indices(moi_model, IOMOI.VariableIndex, IOMOI.GreaterThan{Float64})
    upper_bounds = io_constraint_indices(moi_model, IOMOI.VariableIndex, IOMOI.LessThan{Float64})
    equalities = io_constraint_indices(moi_model, IOMOI.VariableIndex, IOMOI.EqualTo{Float64})

    @test length(zero_one) == 1
    @test IOMOI.get(moi_model, IOMOI.ConstraintFunction(), zero_one[1]) == moi_vars[1]

    @test length(integers) == 3
    @test sort!(collect(IOMOI.get(moi_model, IOMOI.ConstraintFunction(), ci).value for ci in integers)) == [2, 3, 4]

    @test length(intervals) == 1
    @test IOMOI.get(moi_model, IOMOI.ConstraintFunction(), intervals[1]) == moi_vars[1]
    interval = IOMOI.get(moi_model, IOMOI.ConstraintSet(), intervals[1])
    @test interval.lower == 0.0
    @test interval.upper == 1.0

    @test length(upper_bounds) == 1
    @test IOMOI.get(moi_model, IOMOI.ConstraintFunction(), upper_bounds[1]) == moi_vars[2]
    @test IOMOI.get(moi_model, IOMOI.ConstraintSet(), upper_bounds[1]).upper == 5.0

    @test length(lower_bounds) == 1
    @test IOMOI.get(moi_model, IOMOI.ConstraintFunction(), lower_bounds[1]) == moi_vars[3]
    @test IOMOI.get(moi_model, IOMOI.ConstraintSet(), lower_bounds[1]).lower == -2.0

    @test length(equalities) == 1
    @test IOMOI.get(moi_model, IOMOI.ConstraintFunction(), equalities[1]) == moi_vars[4]
    @test IOMOI.get(moi_model, IOMOI.ConstraintSet(), equalities[1]).value == 4.0
end

@testset "build_moi_model exports quadratic constraints and objectives" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 2.0),
        2 => PC.IntVar(0.0, 3.0),
    )
    con = PC.Constraint(
        1,
        PC.QuadExpr(IOQuadTerm[(2.0, 1, 1), (3.0, 1, 2)], IOLinTerm[(5.0, 1)]),
        -1.0,
        9.0,
    )
    obj = PC.QuadExpr(
        IOQuadTerm[(7.0, 2, 2), (11.0, 1, 2)],
        IOLinTerm[(13.0, 2)];
        constant = 17.0,
    )
    model = PC.QPModel(vars, [con], obj, :max)

    moi_model = QIPresolve.build_moi_model(model)

    @test IOMOI.get(moi_model, IOMOI.ObjectiveSense()) == IOMOI.MAX_SENSE
    @test IOMOI.get(moi_model, IOMOI.ObjectiveFunctionType()) == IOMOI.ScalarQuadraticFunction{Float64}

    objective = IOMOI.get(moi_model, IOMOI.ObjectiveFunction{IOMOI.ScalarQuadraticFunction{Float64}}())
    objective_quad = Dict((term.variable_1.value, term.variable_2.value) => term.coefficient for term in objective.quadratic_terms)
    objective_lin = Dict(term.variable.value => term.coefficient for term in objective.affine_terms)
    @test objective.constant == 17.0
    @test objective_quad[(2, 2)] == 14.0
    @test objective_quad[(1, 2)] == 11.0
    @test objective_lin[2] == 13.0

    quad_cons = io_constraint_indices(moi_model, IOMOI.ScalarQuadraticFunction{Float64}, IOMOI.Interval{Float64})
    @test length(quad_cons) == 1
    con_function = IOMOI.get(moi_model, IOMOI.ConstraintFunction(), quad_cons[1])
    con_set = IOMOI.get(moi_model, IOMOI.ConstraintSet(), quad_cons[1])
    con_quad = Dict((term.variable_1.value, term.variable_2.value) => term.coefficient for term in con_function.quadratic_terms)
    con_lin = Dict(term.variable.value => term.coefficient for term in con_function.affine_terms)

    @test con_set.lower == -1.0
    @test con_set.upper == 9.0
    @test con_quad[(1, 1)] == 4.0
    @test con_quad[(1, 2)] == 3.0
    @test con_lin[1] == 5.0
end

@testset "save_moi writes LP interval function constraints as one-sided constraints" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "ranged_quadratic.lp")
        vars = Dict{PC.VarId, PC.IntVar}(
            1 => PC.IntVar(0.0, 2.0),
            2 => PC.IntVar(0.0, 3.0),
        )
        expected_expr = PC.QuadExpr(IOQuadTerm[(2.0, 1, 1), (3.0, 1, 2)], IOLinTerm[(5.0, 1)])
        con = PC.Constraint(1, deepcopy(expected_expr), -1.0, 9.0)
        model = PC.QPModel(vars, [con], io_empty_qe(), :min)

        QIPresolve.ModelIO.save_moi(QIPresolve.build_moi_model(model), lp_path)
        contents = read(lp_path, String)

        @test occursin("_lb:", contents)
        @test occursin("_ub:", contents)
        @test !occursin(r":\s*-?[\d.]+(?:[eE][+-]?\d+)?\s*<=\s*\[", contents)

        roundtrip = io_load_lp_core_model(lp_path)

        @test length(roundtrip.vars) == 2
        @test roundtrip.vars[1] == vars[1]
        @test roundtrip.vars[2] == vars[2]
        @test length(roundtrip.cons) == 2

        lower = only([constraint for constraint in roundtrip.cons if constraint.lhs == -1.0 && constraint.rhs == Inf])
        upper = only([constraint for constraint in roundtrip.cons if constraint.lhs == -Inf && constraint.rhs == 9.0])
        io_assert_same_expr(lower.qe, expected_expr)
        io_assert_same_expr(upper.qe, expected_expr)
    end
end

@testset "build_moi_model handles objective senses and infeasibility" begin
    expected_senses = Dict(
        :min => IOMOI.MIN_SENSE,
        :max => IOMOI.MAX_SENSE,
        :feas => IOMOI.FEASIBILITY_SENSE,
        :undef => IOMOI.FEASIBILITY_SENSE,
    )

    for (sense, moi_sense) in expected_senses
        model = PC.QPModel(Dict{PC.VarId, PC.IntVar}(), PC.Constraint[], io_empty_qe(), sense)
        moi_model = QIPresolve.build_moi_model(model)
        @test IOMOI.get(moi_model, IOMOI.ObjectiveSense()) == moi_sense
    end

    bad_model = PC.QPModel(Dict{PC.VarId, PC.IntVar}(), PC.Constraint[], io_empty_qe(), :invalid)
    @test_throws ArgumentError QIPresolve.build_moi_model(bad_model)

    infeasible_model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(-Inf, Inf)),
        PC.Constraint[],
        io_empty_qe(),
        :min,
    )
    infeasible_model.infeasible = true
    moi_model = QIPresolve.build_moi_model(infeasible_model)
    infeasible_cons = io_constraint_indices(moi_model, IOMOI.ScalarAffineFunction{Float64}, IOMOI.LessThan{Float64})

    @test length(infeasible_cons) == 1
    infeasible_function = IOMOI.get(moi_model, IOMOI.ConstraintFunction(), infeasible_cons[1])
    infeasible_set = IOMOI.get(moi_model, IOMOI.ConstraintSet(), infeasible_cons[1])
    @test isempty(infeasible_function.terms)
    @test infeasible_function.constant == 0.0
    @test infeasible_set.upper == -1.0
end

@testset "build_moi_model selects SCIP backend" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 3.0),
    )
    quadratic_con = PC.Constraint(
        1,
        PC.QuadExpr(IOQuadTerm[(2.0, 1, 2)], IOLinTerm[(1.0, 1)]),
        -Inf,
        4.0,
    )
    linear_obj = PC.QuadExpr(IOQuadTerm[], IOLinTerm[(3.0, 2)]; constant = 5.0)
    linear_model = PC.QPModel(vars, [quadratic_con], linear_obj, :min)

    default_model = QIPresolve.build_moi_model(linear_model)
    @test default_model isa IOMOI.Utilities.Model{Float64}

    scip_model = QIPresolve.build_moi_model(linear_model, :scip)
    @test scip_model isa IOMOI.ModelLike
    @test IOMOI.get(scip_model, IOMOI.ObjectiveSense()) == IOMOI.MIN_SENSE

    quadratic_obj = PC.QuadExpr(IOQuadTerm[(1.0, 1, 1)], IOLinTerm[(3.0, 2)])
    quadratic_model = PC.QPModel(vars, PC.Constraint[], quadratic_obj, :min)
    @test_throws ArgumentError QIPresolve.build_moi_model(quadratic_model, :scip)
    @test_throws ArgumentError QIPresolve.build_moi_model(linear_model, :unknown)
end

@testset "build_moi_model round-trips through from_moi" begin
    affine_vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(-2.0, 5.0),
    )
    affine_con = PC.Constraint(1, PC.QuadExpr(IOQuadTerm[], IOLinTerm[(2.0, 1), (-3.0, 2)]), -4.0, 6.0)
    affine_obj = PC.QuadExpr(IOQuadTerm[], IOLinTerm[(1.5, 1)]; constant = 2.5)
    affine_model = PC.QPModel(affine_vars, [affine_con], affine_obj, :min)

    affine_roundtrip = QIPresolve.ModelIO.build_model(QIPresolve.ModelIO.from_moi(QIPresolve.build_moi_model(affine_model)))
    @test affine_roundtrip.obj_sense == :min
    @test affine_roundtrip.vars[1] == affine_model.vars[1]
    @test affine_roundtrip.vars[2] == affine_model.vars[2]
    @test length(affine_roundtrip.cons) == 1
    io_assert_same_expr(affine_roundtrip.cons[1].qe, affine_model.cons[1].qe)
    @test affine_roundtrip.cons[1].lhs == affine_model.cons[1].lhs
    @test affine_roundtrip.cons[1].rhs == affine_model.cons[1].rhs
    io_assert_same_expr(affine_roundtrip.obj_expr, affine_model.obj_expr)

    quadratic_vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 2.0),
        2 => PC.IntVar(-1.0, 4.0),
    )
    quadratic_con = PC.Constraint(
        1,
        PC.QuadExpr(IOQuadTerm[(2.0, 1, 1), (3.0, 1, 2)], IOLinTerm[(4.0, 2)]),
        -5.0,
        7.0,
    )
    quadratic_obj = PC.QuadExpr(
        IOQuadTerm[(5.0, 2, 2), (6.0, 1, 2)],
        IOLinTerm[(7.0, 1)];
        constant = 8.0,
    )
    quadratic_model = PC.QPModel(quadratic_vars, [quadratic_con], quadratic_obj, :max)

    quadratic_roundtrip = QIPresolve.ModelIO.build_model(QIPresolve.ModelIO.from_moi(QIPresolve.build_moi_model(quadratic_model)))
    @test quadratic_roundtrip.obj_sense == :max
    @test quadratic_roundtrip.vars[1] == quadratic_model.vars[1]
    @test quadratic_roundtrip.vars[2] == quadratic_model.vars[2]
    @test length(quadratic_roundtrip.cons) == 1
    io_assert_same_expr(quadratic_roundtrip.cons[1].qe, quadratic_model.cons[1].qe)
    @test quadratic_roundtrip.cons[1].lhs == quadratic_model.cons[1].lhs
    @test quadratic_roundtrip.cons[1].rhs == quadratic_model.cons[1].rhs
    io_assert_same_expr(quadratic_roundtrip.obj_expr, quadratic_model.obj_expr)
end
