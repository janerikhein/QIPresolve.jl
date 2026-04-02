using Test
using Graphs
using QIPresolve.GraphEmbedding: generate_2_connected_instance, random_2_connected_graph


canonical_edges(g::Graphs.AbstractGraph) = sort(map(e -> (min(src(e), dst(e)), max(src(e), dst(e))), collect(edges(g))))
coord_tuples(coords) = [(p.x, p.y) for p in coords]


@testset "random_2_connected_graph validation" begin
    @test_throws ArgumentError random_2_connected_graph(2)
    @test_throws ArgumentError random_2_connected_graph(3; R = -1)
    @test_throws ArgumentError random_2_connected_graph(5; R = 0)
    @test_throws ArgumentError random_2_connected_graph(5; edge_density = -0.1)
    @test_throws ArgumentError random_2_connected_graph(5; edge_density = 1.1)
end

@testset "random_2_connected_graph edge counts" begin
    g0, coords0 = random_2_connected_graph(6; R = 10, edge_density = 0.0, seed = 11)
    @test nv(g0) == 6
    @test ne(g0) == 6
    @test length(coords0) == 6
    @test isempty(articulation(g0))

    g1, _ = random_2_connected_graph(6; R = 10, edge_density = 1.0, seed = 11)
    @test ne(g1) == 15
    @test isempty(articulation(g1))

    gm, _ = random_2_connected_graph(6; R = 10, edge_density = 0.4, seed = 11)
    @test ne(gm) == round(Int, 6 + 0.4 * (15 - 6))
    @test isempty(articulation(gm))
end

@testset "random_2_connected_graph coordinates and reproducibility" begin
    n = 8
    R = 4
    g1, coords1 = random_2_connected_graph(n; R = R, edge_density = 0.35, seed = 123)
    g2, coords2 = random_2_connected_graph(n; R = R, edge_density = 0.35, seed = 123)

    tuples1 = coord_tuples(coords1)
    tuples2 = coord_tuples(coords2)

    @test tuples1 == tuples2
    @test canonical_edges(g1) == canonical_edges(g2)
    @test allunique(tuples1)
    @test all(t -> -R <= t[1] <= R && -R <= t[2] <= R, tuples1)
    @test isempty(articulation(g1))
end

@testset "generate_2_connected_instance" begin
    model, x, y = generate_2_connected_instance(7; R = 12, edge_density = 0.5, seed = 7)

    @test length(x) == 7
    @test length(y) == 7
    @test model !== nothing
end
