using Test
import QIPresolve.PresolvingCore as PC

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}
const _NEXT_CON_ID = Ref(0)

next_con_id() = (_NEXT_CON_ID[] += 1)

function make_model()
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 2.0),
        2 => PC.IntVar(-1.0, 3.0),
    )

    quad_terms = QuadTerm[(1.0, 1, 1), (0.5, 1, 2)]
    lin_terms = LinTerm[(2.0, 1), (-1.0, 2)]
    con = PC.Constraint(next_con_id(), PC.QuadExpr(quad_terms, lin_terms; constant = 1.2), -3.0, 4.0)

    obj_terms = QuadTerm[(2.0, 1, 2)]
    obj_lin = LinTerm[(1.0, 2)]
    obj = PC.QuadExpr(obj_terms, obj_lin; constant = -0.5)

    return PC.QPModel(vars, [con], obj, :min)
end

@testset "QPModel constructor and variable ops" begin
    model = make_model()
    @test model._max_var_id == 2

    new_id = PC.add_var!(model, PC.IntVar(-2.0, 5.0))
    @test new_id == 3
    @test model.vars[3].lb == -2.0
    @test model.vars[3].ub == 5.0

    PC.set_var_bounds!(model, 3, -1.5, 4.5)
    @test model.vars[3].lb == -1.5
    @test model.vars[3].ub == 4.5
end

@testset "QPModel affine_transform!" begin
    model0 = make_model()
    model = deepcopy(model0)

    PC.affine_transform!(model, 1, 1.3, -0.6)

    con0 = model0.cons[1]
    con = model.cons[1]
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs
    @test isapprox(shift_lhs, shift_rhs; atol = 1.0e-12)
    @test con.qe.constant == 0.0

    for _ in 1:5
        x = randn(2)
        x_sub = copy(x)
        x_sub[1] = 1.3 * x[1] - 0.6

        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1.0e-8)

        obj_before = PC.eval_full(model0.obj_expr, x_sub)
        obj_after = PC.eval_full(model.obj_expr, x)
        @test isapprox(obj_after, obj_before; atol = 1.0e-8)
    end
end

@testset "QPModel lin_transform!" begin
    model0 = make_model()
    model = deepcopy(model0)

    PC.lin_transform!(model, 1, 2, 1.5, -0.7)

    con0 = model0.cons[1]
    con = model.cons[1]
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs
    @test isapprox(shift_lhs, shift_rhs; atol = 1.0e-12)
    @test con.qe.constant == 0.0

    for _ in 1:5
        x = randn(2)
        x_sub = copy(x)
        x_sub[1] = 1.5 * x[1] - 0.7 * x[2]

        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1.0e-8)

        obj_before = PC.eval_full(model0.obj_expr, x_sub)
        obj_after = PC.eval_full(model.obj_expr, x)
        @test isapprox(obj_after, obj_before; atol = 1.0e-8)
    end
end

@testset "QPModel lin_transform! applies each constraint once and supports offsets" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(-2.0, 4.0),
        2 => PC.IntVar(-1.0, 3.0),
        3 => PC.IntVar(0.0, 1.0),
    )
    con1 = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(1.0, 1, 1)], LinTerm[(2.0, 1), (1.5, 2)]),
        -2.0,
        5.0,
    )
    con2 = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(0.5, 1, 2)], LinTerm[(1.0, 1), (-2.0, 2)]),
        -3.0,
        4.0,
    )
    obj = PC.QuadExpr(QuadTerm[(1.0, 1, 2)], LinTerm[(3.0, 1), (1.0, 2)])
    model0 = PC.QPModel(vars, [con1, con2], obj, :min)
    model = deepcopy(model0)

    PC.lin_transform!(model, 1, 3, 2.0, -1.0; c = 1.0)

    for con_idx in eachindex(model.cons)
        con0 = model0.cons[con_idx]
        con = model.cons[con_idx]
        shift_lhs = con0.lhs - con.lhs
        shift_rhs = con0.rhs - con.rhs
        @test isapprox(shift_lhs, shift_rhs; atol = 1.0e-12)
        @test con.qe.constant == 0.0
    end

    @test 3 in collect(PC.vars(model.cons[1]))
    @test 3 in collect(PC.vars(model.cons[2]))
    @test 3 in collect(PC.vars(model.obj_expr))

    for _ in 1:5
        x = randn(3)
        x_sub = copy(x)
        x_sub[1] = 2.0 * x[1] - x[3] + 1.0

        for con_idx in eachindex(model.cons)
            con0 = model0.cons[con_idx]
            con = model.cons[con_idx]
            shift = con0.lhs - con.lhs
            val_before = PC.eval_full(con0.qe, x_sub)
            val_after = PC.eval_full(con.qe, x)
            @test isapprox(val_after, val_before - shift; atol = 1.0e-8)
        end

        obj_before = PC.eval_full(model0.obj_expr, x_sub)
        obj_after = PC.eval_full(model.obj_expr, x)
        @test isapprox(obj_after, obj_before; atol = 1.0e-8)
    end
end

@testset "QPModel var_bound_shift!" begin
    model0 = make_model()
    model = deepcopy(model0)

    shift = 0.5
    PC.var_bound_shift!(model, 1, shift)

    @test model.vars[1].lb == model0.vars[1].lb - shift
    @test model.vars[1].ub == model0.vars[1].ub - shift

    con0 = model0.cons[1]
    con = model.cons[1]
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs
    @test isapprox(shift_lhs, shift_rhs; atol = 1.0e-12)
    @test con.qe.constant == 0.0

    for _ in 1:5
        x = randn(2)
        x_sub = copy(x)
        x_sub[1] = x[1] + shift

        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1.0e-8)

        obj_before = PC.eval_full(model0.obj_expr, x_sub)
        obj_after = PC.eval_full(model.obj_expr, x)
        @test isapprox(obj_after, obj_before; atol = 1.0e-8)
    end
end

@testset "QPModel scale_constraints_gcd!" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 3.0),
        2 => PC.IntVar(-1.0, 4.0),
    )
    con1 = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(6.0, 1, 1), (8.0, 1, 2)], LinTerm[(4.0, 1), (-2.0, 2)]),
        -4.0,
        10.0,
    )
    con2 = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(4.0, 2, 2)], LinTerm[(1.5, 2)]; constant = 1.0),
        -1.0,
        5.0,
    )
    obj = PC.QuadExpr(QuadTerm[], LinTerm[])
    model = PC.QPModel(vars, [con1, con2], obj, :min)

    @test PC.scale_constraints_gcd!(model) == 1
    @test !model.infeasible

    @test model.cons[1].lhs == -2.0
    @test model.cons[1].rhs == 5.0
    @test PC.get_quad_coeff(model.cons[1].qe, 1, 1) == 3.0
    @test PC.get_quad_coeff(model.cons[1].qe, 1, 2) == 4.0
    @test PC.get_lin_coeff(model.cons[1].qe, 1) == 2.0
    @test PC.get_lin_coeff(model.cons[1].qe, 2) == -1.0

    @test model.cons[2].qe.constant == 0.0
    @test model.cons[2].lhs == -2.0
    @test model.cons[2].rhs == 4.0
    @test PC.get_quad_coeff(model.cons[2].qe, 2, 2) == 4.0
    @test PC.get_lin_coeff(model.cons[2].qe, 2) == 1.5
end

@testset "QPModel scale_constraints_gcd! detects infeasibility" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 2.0))
    con = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(2.0, 1, 1)], LinTerm[]),
        1.0,
        1.0,
    )
    model = PC.QPModel(vars, [con], PC.QuadExpr(QuadTerm[], LinTerm[]), :min)

    @test PC.scale_constraints_gcd!(model) == 1
    @test model.infeasible
    @test model.cons[1].lhs == 1.0
    @test model.cons[1].rhs == 0.0
end

@testset "QPModel scale_constraints_gcd! respects infeasible flag" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 2.0))
    con = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(4.0, 1, 1)], LinTerm[(2.0, 1)]),
        0.0,
        6.0,
    )
    model = PC.QPModel(vars, [con], PC.QuadExpr(QuadTerm[], LinTerm[]), :min)
    model.infeasible = true

    @test PC.scale_constraints_gcd!(model) == 0
    @test PC.get_quad_coeff(model.cons[1].qe, 1, 1) == 4.0
    @test PC.get_lin_coeff(model.cons[1].qe, 1) == 2.0
end

@testset "QPModel fix_vars! handles singleton linear constraints with negative coefficients" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0))
    ranged = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[], LinTerm[(-3.0, 1)]),
        -4.0,
        0.0,
    )
    ranged_model = PC.QPModel(deepcopy(vars), [ranged], PC.QuadExpr(QuadTerm[], LinTerm[]), :min)

    PC.fix_vars!(ranged_model)

    @test !ranged_model.infeasible
    @test ranged_model.vars[1] == PC.IntVar(0.0, 1.0)
    @test isempty(ranged_model.cons)

    equality = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[], LinTerm[(-2.0, 1)]),
        -2.0,
        -2.0,
    )
    equality_model = PC.QPModel(deepcopy(vars), [equality], PC.QuadExpr(QuadTerm[], LinTerm[]), :min)

    PC.fix_vars!(equality_model)

    @test !equality_model.infeasible
    @test equality_model.vars[1] == PC.IntVar(1.0, 1.0)
    @test isempty(equality_model.cons)

    infeasible = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[], LinTerm[(-2.0, 1)]),
        -5.0,
        -5.0,
    )
    infeasible_model = PC.QPModel(deepcopy(vars), [infeasible], PC.QuadExpr(QuadTerm[], LinTerm[]), :min)

    PC.fix_vars!(infeasible_model)

    @test infeasible_model.infeasible
end
