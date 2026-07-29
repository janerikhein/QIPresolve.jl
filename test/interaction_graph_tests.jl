using Test
using Graphs: has_edge, ne, nv
import QIPresolve.PresolvingCore as PC

const IGQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const IGLinTerm = Tuple{Float64, PC.VarId}

@testset "InteractionGraph builds constraint bilinear residue graph" begin
    con = PC.Constraint(
        1,
        PC.QuadExpr(
            IGQuadTerm[
                (4.0, 10, 2),
                (-5.0, 2, 7),
                (6.0, 7, 10),
                (9.0, 2, 2),
            ],
            IGLinTerm[(3.0, 5)],
        ),
        -10.0,
        10.0,
    )

    ig = PC.InteractionGraph(con, 3)

    @test ig.modulus == 3
    @test ig.pos_to_var_id == [2, 5, 7, 10]
    @test ig.var_id_to_pos == Dict(2 => 1, 5 => 2, 7 => 3, 10 => 4)
    @test nv(ig.graph) == 4
    @test ne(ig.graph) == 2
    @test ig.lin_coeffs == Dict{PC.VarId, Int}()
    @test ig.quad_coeffs == Dict{Tuple{PC.VarId, PC.VarId}, Int}((2, 10) => 1, (2, 7) => 1)
    @test has_edge(ig.graph, ig.var_id_to_pos[2], ig.var_id_to_pos[10])
    @test has_edge(ig.graph, ig.var_id_to_pos[2], ig.var_id_to_pos[7])
    @test !has_edge(ig.graph, ig.var_id_to_pos[7], ig.var_id_to_pos[10])
    @test !has_edge(ig.graph, ig.var_id_to_pos[2], ig.var_id_to_pos[2])
    @test !has_edge(ig.graph, ig.var_id_to_pos[5], ig.var_id_to_pos[7])
end

@testset "InteractionGraph handles modulus one and validates inputs" begin
    con = PC.Constraint(
        2,
        PC.QuadExpr(
            IGQuadTerm[(1.0, 1, 2), (-2.0, 2, 3), (3.0, 1, 3)],
            IGLinTerm[],
        ),
        0.0,
        1.0,
    )

    ig = PC.InteractionGraph(con, 1)

    @test ig.pos_to_var_id == [1, 2, 3]
    @test nv(ig.graph) == 3
    @test ne(ig.graph) == 0

    @test_throws ArgumentError PC.InteractionGraph(con, 0)
    @test_throws ArgumentError PC.InteractionGraph(con, -2)
end

@testset "InteractionGraph rejects non-integer bilinear coefficients" begin
    con = PC.Constraint(
        3,
        PC.QuadExpr(IGQuadTerm[(1.5, 1, 2)], IGLinTerm[]),
        0.0,
        1.0,
    )

    @test_throws ArgumentError PC.InteractionGraph(con, 2)
end

@testset "InteractionGraph rejects non-integer singleton coefficients" begin
    lin_con = PC.Constraint(
        4,
        PC.QuadExpr(IGQuadTerm[], IGLinTerm[(1.5, 1)]),
        0.0,
        1.0,
    )
    diag_con = PC.Constraint(
        5,
        PC.QuadExpr(IGQuadTerm[(1.5, 1, 1)], IGLinTerm[]),
        0.0,
        1.0,
    )

    @test_throws ArgumentError PC.InteractionGraph(lin_con, 2)
    @test_throws ArgumentError PC.InteractionGraph(diag_con, 2)
end

@testset "decompose returns residual interaction components" begin
    con = PC.Constraint(
        6,
        PC.QuadExpr(
            IGQuadTerm[
                (7.0, 8, 1),
                (-4.0, 3, 8),
                (10.0, 1, 3),
                (6.0, 1, 1),
                (-3.0, 6, 6),
            ],
            IGLinTerm[(12.0, 1), (-2.0, 8), (7.0, 4), (10.0, 6), (5.0, 7)],
        ),
        -20.0,
        20.0,
    )

    ig = PC.InteractionGraph(con, 5)

    @test ig.pos_to_var_id == [1, 3, 4, 6, 7, 8]
    @test ig.lin_coeffs == Dict{PC.VarId, Int}(1 => 2, 8 => 3, 4 => 2)
    @test ig.quad_coeffs == Dict{Tuple{PC.VarId, PC.VarId}, Int}(
        (1, 8) => 2,
        (3, 8) => 1,
        (1, 1) => 1,
        (6, 6) => 2,
    )

    components = PC.decompose(ig)

    @test length(components) == 3

    non_singleton = components[1]
    @test non_singleton isa PC.NonSingleton
    @test non_singleton.pos_to_var_id == [1, 3, 8]
    @test non_singleton.var_id_to_pos == Dict(1 => 1, 3 => 2, 8 => 3)
    @test nv(non_singleton.graph) == 3
    @test ne(non_singleton.graph) == 2
    @test has_edge(non_singleton.graph, 1, 3)
    @test has_edge(non_singleton.graph, 2, 3)
    @test non_singleton.lin_coeffs == Dict{PC.VarId, Int}(1 => 2, 8 => 3)
    @test non_singleton.quad_coeffs == Dict{Tuple{PC.VarId, PC.VarId}, Int}(
        (1, 8) => 2,
        (3, 8) => 1,
        (1, 1) => 1,
    )

    lin_singleton = components[2]
    @test lin_singleton isa PC.LinSingleton
    @test lin_singleton.var_id == 4
    @test lin_singleton.lin_coeff == 2

    quad_singleton = components[3]
    @test quad_singleton isa PC.QuadSingleton
    @test quad_singleton.var_id == 6
    @test quad_singleton.lin_coeff == 0
    @test quad_singleton.quad_coeff == 2
end
