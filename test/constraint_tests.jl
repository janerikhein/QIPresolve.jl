using Test
import QIPresolve.PresolvingCore as PC

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}

@testset "Constraint normalize!" begin
    quad_terms = QuadTerm[(2.0, 1, 1), (1.0, 1, 2)]
    lin_terms = LinTerm[(3.0, 1), (-1.0, 2)]
    qe = PC.QuadExpr(quad_terms, lin_terms; constant = 4.5)
    con = PC.Constraint(qe, -2.0, 7.0)

    x = randn(2)
    lhs_before = con.lhs
    rhs_before = con.rhs
    eval_before = PC.eval_full(con.qe, x)

    PC.normalize!(con)

    @test con.qe.constant == 0.0
    @test con.lhs == lhs_before - 4.5
    @test con.rhs == rhs_before - 4.5
    @test isapprox(PC.eval_full(con.qe, x), eval_before - 4.5; atol = 1e-12)
end

@testset "Constraint symmetrize!" begin
    quad_terms = QuadTerm[(1.0, 1, 2), (3.0, 2, 1), (2.0, 1, 1)]
    lin_terms = LinTerm[(4.0, 1), (-2.0, 2)]
    qe = PC.QuadExpr(quad_terms, lin_terms; constant = 1.0)
    con = PC.Constraint(qe, -3.0, 5.0)

    x = randn(2)
    eval_before = PC.eval_full(con.qe, x)
    lhs_before = con.lhs
    rhs_before = con.rhs
    lin_before = copy(PC.lin(con.qe))
    quad_before = copy(PC.quad(con.qe))

    PC.symmetrize!(con)

    @test con.lhs == 2 * lhs_before
    @test con.rhs == 2 * rhs_before
    @test con.qe.constant == 2.0
    @test isapprox(PC.lin(con.qe), 2 * lin_before; atol = 1e-12)
    @test isapprox(PC.quad(con.qe), quad_before + transpose(quad_before); atol = 1e-12)
    @test isapprox(PC.eval_full(con.qe, x), 2 * eval_before; atol = 1e-12)
end

@testset "Constraint affine_transform!" begin
    quad_terms = QuadTerm[(1.0, 1, 1), (2.0, 1, 2), (3.0, 2, 1), (4.0, 2, 2)]
    lin_terms = LinTerm[(5.0, 1), (-2.0, 2)]
    con0 = PC.Constraint(PC.QuadExpr(quad_terms, lin_terms; constant = 1.5), -3.0, 7.0)

    con = deepcopy(con0)
    PC.affine_transform!(con, 1, 1.2, -0.7)

    @test con.qe.constant == 0.0
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs
    @test isapprox(shift_lhs, shift_rhs; atol = 1e-12)

    for _ in 1:5
        x = randn(2)
        x_sub = copy(x)
        x_sub[1] = 1.2 * x[1] - 0.7
        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1e-8)
    end
end

@testset "Constraint lin_transform!" begin
    quad_terms = QuadTerm[
        (1.0, 1, 1),
        (2.0, 1, 2),
        (3.0, 2, 1),
        (4.0, 2, 2),
        (5.0, 1, 3),
        (6.0, 3, 1),
        (7.0, 2, 3),
        (8.0, 3, 2),
        (9.0, 3, 3),
    ]
    lin_terms = LinTerm[(1.0, 1), (2.0, 2), (3.0, 3)]
    con0 = PC.Constraint(PC.QuadExpr(quad_terms, lin_terms; constant = -2.5), -4.0, 6.0)

    con = deepcopy(con0)
    PC.lin_transform!(con, 1, 2, 1.5, -0.7)

    @test con.qe.constant == 0.0
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs
    @test isapprox(shift_lhs, shift_rhs; atol = 1e-12)

    for _ in 1:5
        x = randn(3)
        x_sub = copy(x)
        x_sub[1] = 1.5 * x[1] - 0.7 * x[2]
        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1e-8)
    end
end
