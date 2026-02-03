using Random
using Test
import QIPresolve.PresolvingCore as PC

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}

empty_qe() = PC.QuadExpr(QuadTerm[], LinTerm[])

function assert_same_expr(a::PC.QuadExpr, b::PC.QuadExpr; atol::Float64 = 1e-12)
    @test a.nvars == b.nvars
    @test a.cap == b.cap
    @test a.pos_to_var[1:a.nvars] == b.pos_to_var[1:b.nvars]
    @test a.var_to_pos == b.var_to_pos
    @test a.perm == b.perm
    @test isapprox(a.constant, b.constant; atol = atol)
    @test isapprox(a.lin_buf, b.lin_buf; atol = atol)
    @test isapprox(a.quad_buf, b.quad_buf; atol = atol)
end

@testset "QuadExpr constructor and views" begin
    quad_terms = QuadTerm[(1.5, 2, 1), (2.0, 2, 1), (-0.5, 1, 1)]
    lin_terms = LinTerm[(2.0, 1), (-1.0, 2), (3.0, 3)]
    qe = PC.QuadExpr(quad_terms, lin_terms; constant = 1.2, buf_size_factor = 1.5)

    @test qe.nvars == 3
    @test qe.constant == 1.2

    @test qe.var_to_pos[2] == 1
    @test qe.var_to_pos[1] == 2
    @test qe.var_to_pos[3] == 3

    @test size(PC.quad(qe)) == (3, 3)
    @test size(PC.lin(qe)) == (3,)
    @test collect(PC.vars(qe)) == [2, 1, 3]

    p2 = qe.var_to_pos[2]
    p1 = qe.var_to_pos[1]
    p3 = qe.var_to_pos[3]

    @test PC.quad(qe)[p2, p1] == 3.5
    @test PC.quad(qe)[p1, p1] == -0.5

    @test PC.lin(qe)[p1] == 2.0
    @test PC.lin(qe)[p2] == -1.0
    @test PC.lin(qe)[p3] == 3.0

    @test_throws AssertionError PC.QuadExpr(QuadTerm[], LinTerm[]; buf_size_factor = 0.5)
end

@testset "QuadExpr coefficient updates" begin
    qe = empty_qe()

    @test PC.add_var!(qe, 10) == 1
    @test PC.add_var!(qe, 11) == 2
    @test PC.add_var!(qe, 10) == 1

    @test PC.set_lin_coeff!(qe, 12, 1.0) == false
    @test PC.set_lin_coeff!(qe, 10, 2.5) == true
    @test PC.add_lin_coeff!(qe, 10, 0.5) == true
    @test PC.add_lin_coeff!(qe, 12, 1.0) == false

    p10 = qe.var_to_pos[10]
    p11 = qe.var_to_pos[11]

    @test PC.lin(qe)[p10] == 3.0
    @test PC.get_lin_coeff(qe, 10) == 3.0
    @test PC.get_lin_coeff(qe, 12) == 0.0

    @test PC.set_quad_coeff!(qe, 10, 11, 4.0) == true
    @test PC.quad(qe)[p10, p11] == 4.0
    @test PC.quad(qe)[p11, p10] == 0.0
    @test PC.get_quad_coeff(qe, 10, 11) == 4.0
    @test PC.get_quad_coeff(qe, 11, 10) == 0.0
    @test PC.get_quad_coeff(qe, 10, 12) == 0.0

    @test PC.add_quad_coeff!(qe, 10, 11, 1.5; sym = true) == true
    @test PC.quad(qe)[p10, p11] == 5.5
    @test PC.quad(qe)[p11, p10] == 1.5
    @test PC.get_quad_coeff(qe, 10, 11) == 5.5
    @test PC.get_quad_coeff(qe, 11, 10) == 1.5

    @test PC.set_quad_coeff!(qe, 10, 11, 7.0; sym = true) == true
    @test PC.quad(qe)[p10, p11] == 7.0
    @test PC.quad(qe)[p11, p10] == 7.0

    @test PC.set_quad_coeff!(qe, 10, 12, 1.0) == false
    @test PC.add_quad_coeff!(qe, 10, 12, 1.0) == false
end

@testset "QuadExpr add/remove and perm mapping" begin
    qe = empty_qe()
    PC.add_var!(qe, 1)
    PC.add_var!(qe, 2)
    PC.add_var!(qe, 3)

    PC.set_lin_coeff!(qe, 1, 1.0)
    PC.set_lin_coeff!(qe, 2, 2.0)
    PC.set_lin_coeff!(qe, 3, 3.0)

    PC.set_quad_coeff!(qe, 1, 3, 30.0)
    PC.set_quad_coeff!(qe, 1, 2, 10.0)
    PC.set_quad_coeff!(qe, 2, 3, 20.0)

    pos2 = qe.var_to_pos[2]
    phys2 = qe.perm[pos2]

    @test PC.remove_var!(qe, 2) == true
    @test qe.nvars == 2
    @test !haskey(qe.var_to_pos, 2)
    @test sort(qe.pos_to_var[1:qe.nvars]) == [1, 3]
    @test collect(PC.vars(qe)) == [1, 3]

    p1 = qe.var_to_pos[1]
    p3 = qe.var_to_pos[3]

    @test PC.lin(qe)[p1] == 1.0
    @test PC.lin(qe)[p3] == 3.0
    @test PC.quad(qe)[p1, p3] == 30.0
    @test PC.quad(qe)[p3, p1] == 0.0
    @test PC.get_lin_coeff(qe, 2) == 0.0
    @test PC.get_quad_coeff(qe, 2, 1) == 0.0

    freed_phys = qe.perm[qe.nvars + 1]
    @test freed_phys == phys2
    @test qe.lin_buf[freed_phys] == 0.0
    @test all(iszero, qe.quad_buf[freed_phys, :])
    @test all(iszero, qe.quad_buf[:, freed_phys])

    PC.add_var!(qe, 4)
    @test qe.var_to_pos[4] == 3
    p4 = qe.var_to_pos[4]
    @test PC.lin(qe)[p4] == 0.0
    @test PC.quad(qe)[p4, p1] == 0.0
    @test collect(PC.vars(qe)) == [1, 3, 4]
end

@testset "QuadExpr capacity growth" begin
    qe = empty_qe()
    @test qe.cap == 1

    PC.add_var!(qe, 1)
    @test qe.cap == 1

    PC.add_var!(qe, 2)
    @test qe.cap == 2

    PC.add_var!(qe, 3)
    @test qe.cap == 4
end

@testset "QuadExpr affine_transform! and invert_affine" begin
    @test PC.invert_affine(2.0, -3.0) == (0.5, 1.5)
    @test_throws ArgumentError PC.invert_affine(0.0, 1.0)

    quad_terms = QuadTerm[(1.0, 1, 1), (2.0, 1, 2), (3.0, 2, 1), (4.0, 2, 2)]
    lin_terms = LinTerm[(5.0, 1), (6.0, 2)]
    qe0 = PC.QuadExpr(quad_terms, lin_terms; constant = 7.0, buf_size_factor = 2.0)

    qe_missing = deepcopy(qe0)
    PC.affine_transform!(qe_missing, 99, 1.2, -0.4)
    assert_same_expr(qe_missing, qe0)

    qe_identity = deepcopy(qe0)
    PC.affine_transform!(qe_identity, 1, 1.0, 0.0)
    assert_same_expr(qe_identity, qe0)

    qe_scale0 = deepcopy(qe0)
    PC.affine_transform!(qe_scale0, 1, 0.0, 2.0)
    @test isapprox(qe_scale0.constant, 7.0 + 2.0 * 5.0 + 4.0 * 1.0; atol = 1e-12)
    @test isapprox(PC.get_lin_coeff(qe_scale0, 1), 0.0; atol = 1e-12)
    @test isapprox(PC.get_lin_coeff(qe_scale0, 2), 6.0 + 2.0 * (2.0 + 3.0); atol = 1e-12)
    @test isapprox(PC.get_quad_coeff(qe_scale0, 1, 1), 0.0; atol = 1e-12)
    @test isapprox(PC.get_quad_coeff(qe_scale0, 1, 2), 0.0; atol = 1e-12)
    @test isapprox(PC.get_quad_coeff(qe_scale0, 2, 1), 0.0; atol = 1e-12)
    @test isapprox(PC.get_quad_coeff(qe_scale0, 2, 2), 4.0; atol = 1e-12)

    qe_invert = deepcopy(qe0)
    PC.affine_transform!(qe_invert, 2, 1.7, -0.6)
    PC.affine_transform!(qe_invert, 2, 1.7, -0.6; invert = true)
    assert_same_expr(qe_invert, qe0)
end

@testset "QuadExpr lin_transform! and invert_lin" begin
    @test PC.invert_lin(2.0, -3.0) == (0.5, 1.5)
    @test_throws ArgumentError PC.invert_lin(0.0, 1.0)

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
    qe0 = PC.QuadExpr(quad_terms, lin_terms; constant = -2.5, buf_size_factor = 2.0)

    qe_missing = deepcopy(qe0)
    PC.lin_transform!(qe_missing, 99, 1, 1.2, -0.4)
    assert_same_expr(qe_missing, qe0)

    qe_missing_other = deepcopy(qe0)
    PC.lin_transform!(qe_missing_other, 1, 99, 1.2, -0.4)
    assert_same_expr(qe_missing_other, qe0)

    qe_same = deepcopy(qe0)
    PC.lin_transform!(qe_same, 1, 1, 1.2, -0.4)
    assert_same_expr(qe_same, qe0)

    for (a, b) in ((1.5, -0.7), (0.0, 2.0))
        qe = deepcopy(qe0)
        PC.lin_transform!(qe, 1, 2, a, b)
        for _ in 1:5
            x = randn(3)
            x_sub = copy(x)
            x_sub[1] = a * x[1] + b * x[2]
            @test isapprox(PC.eval_full(qe, x), PC.eval_full(qe0, x_sub); atol = 1e-8)
        end
    end

    qe_inv = deepcopy(qe0)
    PC.lin_transform!(qe_inv, 1, 2, 1.5, -0.7)
    PC.lin_transform!(qe_inv, 1, 2, 1.5, -0.7; invert = true)
    assert_same_expr(qe_inv, qe0)
end
