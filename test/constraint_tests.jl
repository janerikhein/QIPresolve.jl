using Test
import QIPresolve.PresolvingCore as PC

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}
const _NEXT_CON_ID = Ref(0)

next_con_id() = (_NEXT_CON_ID[] += 1)

@testset "Constraint normalize!" begin
    quad_terms = QuadTerm[(2.0, 1, 1), (1.0, 1, 2)]
    lin_terms = LinTerm[(3.0, 1), (-1.0, 2)]
    qe = PC.QuadExpr(quad_terms, lin_terms; constant = 4.5)
    con = PC.Constraint(next_con_id(), qe, -2.0, 7.0)

    x = randn(2)
    lhs_before = con.lhs
    rhs_before = con.rhs
    eval_before = PC.eval_full(con.qe, x)

    PC.normalize!(con)

    @test con.qe.constant == 0.0
    @test con.lhs == lhs_before - 4.5
    @test con.rhs == rhs_before - 4.5
    @test isapprox(PC.eval_full(con.qe, x), eval_before - 4.5; atol = 1.0e-12)
end

@testset "Constraint symmetrize!" begin
    quad_terms = QuadTerm[(1.0, 1, 2), (3.0, 2, 1), (2.0, 1, 1)]
    lin_terms = LinTerm[(4.0, 1), (-2.0, 2)]
    qe = PC.QuadExpr(quad_terms, lin_terms; constant = 1.0)
    con = PC.Constraint(next_con_id(), qe, -3.0, 5.0)

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
    @test isapprox(PC.lin(con.qe), 2 * lin_before; atol = 1.0e-12)
    @test isapprox(PC.quad(con.qe), quad_before + transpose(quad_before); atol = 1.0e-12)
    @test isapprox(PC.eval_full(con.qe, x), 2 * eval_before; atol = 1.0e-12)
end

@testset "Constraint affine_transform!" begin
    quad_terms = QuadTerm[(1.0, 1, 1), (2.0, 1, 2), (3.0, 2, 1), (4.0, 2, 2)]
    lin_terms = LinTerm[(5.0, 1), (-2.0, 2)]
    con0 = PC.Constraint(next_con_id(), PC.QuadExpr(quad_terms, lin_terms; constant = 1.5), -3.0, 7.0)

    con = deepcopy(con0)
    PC.affine_transform!(con, 1, 1.2, -0.7)

    @test con.qe.constant == 0.0
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs
    @test isapprox(shift_lhs, shift_rhs; atol = 1.0e-12)

    for _ in 1:5
        x = randn(2)
        x_sub = copy(x)
        x_sub[1] = 1.2 * x[1] - 0.7
        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1.0e-8)
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
    con0 = PC.Constraint(next_con_id(), PC.QuadExpr(quad_terms, lin_terms; constant = -2.5), -4.0, 6.0)

    con = deepcopy(con0)
    PC.lin_transform!(con, 1, 2, 1.5, -0.7)

    @test con.qe.constant == 0.0
    shift_lhs = con0.lhs - con.lhs
    shift_rhs = con0.rhs - con.rhs
    @test isapprox(shift_lhs, shift_rhs; atol = 1.0e-12)

    for _ in 1:5
        x = randn(3)
        x_sub = copy(x)
        x_sub[1] = 1.5 * x[1] - 0.7 * x[2]
        val_before = PC.eval_full(con0.qe, x_sub)
        val_after = PC.eval_full(con.qe, x)
        @test isapprox(val_after, val_before - shift_lhs; atol = 1.0e-8)
    end
end

@testset "Constraint is_integer" begin
    integral_con = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(4.0, 1, 1), (6.0, 1, 2)], LinTerm[(2.0, 1), (-3.0, 2)]),
        -5.0,
        7.0,
    )
    @test PC.is_integer(integral_con)

    fractional_linear = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(4.0, 1, 1)], LinTerm[(1.5, 1)]),
        -2.0,
        3.0,
    )
    @test !PC.is_integer(fractional_linear)

    fractional_diagonal = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(2.5, 1, 1)], LinTerm[]),
        -2.0,
        3.0,
    )
    @test !PC.is_integer(fractional_diagonal)

    odd_bilinear = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(3.0, 1, 2)], LinTerm[]),
        -2.0,
        3.0,
    )
    @test !PC.is_integer(odd_bilinear)

    fractional_normalized_bound = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[], LinTerm[(1.0, 1)]; constant = 0.5),
        1.0,
        3.0,
    )
    @test !PC.is_integer(fractional_normalized_bound)

    normalized_integer_bounds = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[], LinTerm[(1.0, 1)]; constant = 0.5),
        0.5,
        1.5,
    )
    @test PC.is_integer(normalized_integer_bounds)

    infinite_bound = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(2.0, 1, 1)], LinTerm[(1.0, 1)]),
        -Inf,
        4.0,
    )
    @test PC.is_integer(infinite_bound)
end

@testset "Constraint scale_gcd!" begin
    con = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(
            QuadTerm[(6.0, 1, 1), (8.0, 1, 2), (10.0, 2, 2)],
            LinTerm[(4.0, 1), (-2.0, 2)];
            constant = 2.0,
        ),
        4.0,
        18.0,
    )

    @test PC.scale_gcd!(con)
    @test con.qe.constant == 0.0
    @test con.lhs == 1.0
    @test con.rhs == 8.0
    @test PC.get_quad_coeff(con.qe, 1, 1) == 3.0
    @test PC.get_quad_coeff(con.qe, 1, 2) == 4.0
    @test PC.get_quad_coeff(con.qe, 2, 2) == 5.0
    @test PC.get_lin_coeff(con.qe, 1) == 2.0
    @test PC.get_lin_coeff(con.qe, 2) == -1.0

    fractional = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(4.0, 1, 1)], LinTerm[(1.5, 1)]; constant = 1.0),
        -2.0,
        3.0,
    )
    @test !PC.scale_gcd!(fractional)
    @test fractional.qe.constant == 0.0
    @test fractional.lhs == -3.0
    @test fractional.rhs == 2.0
    @test PC.get_quad_coeff(fractional.qe, 1, 1) == 4.0
    @test PC.get_lin_coeff(fractional.qe, 1) == 1.5

    unit_gcd = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(3.0, 1, 1), (2.0, 1, 2)], LinTerm[(5.0, 1)]),
        -1.0,
        4.0,
    )
    @test !PC.scale_gcd!(unit_gcd)
    @test PC.get_quad_coeff(unit_gcd.qe, 1, 1) == 3.0
    @test PC.get_quad_coeff(unit_gcd.qe, 1, 2) == 2.0
    @test PC.get_lin_coeff(unit_gcd.qe, 1) == 5.0

    infinite_bound = PC.Constraint(
        next_con_id(),
        PC.QuadExpr(QuadTerm[(4.0, 1, 1)], LinTerm[(2.0, 1)]),
        -Inf,
        6.0,
    )
    @test PC.scale_gcd!(infinite_bound)
    @test infinite_bound.lhs == -Inf
    @test infinite_bound.rhs == 3.0
end
