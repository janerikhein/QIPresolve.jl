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
