using Random
using Test
using QIPresolve

const QE = QIPresolve

function rand_expr(rng::AbstractRNG, nvars::Int, nlin::Int, nquad::Int)
    e = QE.QuadExpr()
    for _ in 1:nlin
        i = rand(rng, 1:nvars)
        QE.add_lin!(e, i, randn(rng))
    end
    for _ in 1:nquad
        i = rand(rng, 1:nvars)
        j = rand(rng, 1:nvars)
        QE.add_quad!(e, i, j, randn(rng))
    end
    QE.add_constant!(e, randn(rng))
    return e
end

@testset "quadexpr canonicalization" begin
    e = QE.QuadExpr()
    QE.add_quad!(e, 3, 1, 2.0)
    @test haskey(e.quad, (1, 3))
    @test !haskey(e.quad, (3, 1))
end

@testset "quadexpr algebra" begin
    rng = MersenneTwister(0)
    for _ in 1:50
        nvars = 6
        e = rand_expr(rng, nvars, 6, 8)
        f = rand_expr(rng, nvars, 5, 7)
        x = randn(rng, nvars)
        alpha = randn(rng)

        e_scaled = copy(e)
        QE.scale!(e_scaled, alpha)
        @test isapprox(QE.eval_expr(e_scaled, x), alpha * QE.eval_expr(e, x); atol = 1.0e-8)

        e_sum = e + f
        @test isapprox(QE.eval_expr(e_sum, x), QE.eval_expr(e, x) + QE.eval_expr(f, x); atol = 1.0e-8)
    end
end

@testset "quadexpr fix_var" begin
    rng = MersenneTwister(1)
    for _ in 1:50
        nvars = 6
        e = rand_expr(rng, nvars, 6, 8)
        i = rand(rng, 1:nvars)
        v = randn(rng)
        x = randn(rng, nvars)
        x2 = copy(x)
        x2[i] = v

        lhs = QE.eval_expr(e, x2)
        e_fix = QE.fix_var(e, i, v)
        rhs = QE.eval_expr(e_fix, x)
        @test isapprox(lhs, rhs; atol = 1.0e-8)
    end
end

@testset "quadexpr substitute_affine" begin
    rng = MersenneTwister(2)
    for _ in 1:50
        nvars = 7
        e = rand_expr(rng, nvars, 6, 8)
        old = rand(rng, 1:nvars)
        new = rand(rng, 1:nvars)
        while new == old
            new = rand(rng, 1:nvars)
        end
        alpha = randn(rng)
        beta = randn(rng)
        x = randn(rng, nvars)
        x2 = copy(x)
        x2[old] = alpha * x2[new] + beta

        lhs = QE.eval_expr(e, x2)
        e_sub = QE.substitute_affine(e, old, alpha, new, beta)
        rhs = QE.eval_expr(e_sub, x)
        @test isapprox(lhs, rhs; atol = 1.0e-8)
    end
end

@testset "quadexpr substitute_affine alpha zero" begin
    rng = MersenneTwister(3)
    for _ in 1:50
        nvars = 6
        e = rand_expr(rng, nvars, 6, 8)
        old = rand(rng, 1:nvars)
        new = rand(rng, 1:nvars)
        while new == old
            new = rand(rng, 1:nvars)
        end
        beta = randn(rng)
        x = randn(rng, nvars)

        e_sub = QE.substitute_affine(e, old, 0.0, new, beta)
        e_fix = QE.fix_var(e, old, beta)
        @test isapprox(QE.eval_expr(e_sub, x), QE.eval_expr(e_fix, x); atol = 1.0e-8)
    end
end
