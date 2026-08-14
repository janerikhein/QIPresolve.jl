using Test
using Graphs
using JuMP: backend
import JuMP
import MathOptInterface as MOI
import QIPresolve
import QIPresolve.PresolvingCore as PC
import QIPresolve.InstanceGeneration as IG
using QIPresolve.InstanceGeneration:
    generate_2_connected_instance,
    generate_globally_rigid_instance,
    generate_laman_instance,
    generate_random_qip_model,
    random_2_connected_graph,
    random_globally_rigid_graph,
    random_laman_graph


canonical_edges(g::Graphs.AbstractGraph) = sort(map(e -> (min(src(e), dst(e)), max(src(e), dst(e))), collect(edges(g))))
coord_tuples(coords) = [(p.x, p.y) for p in coords]
random_qip_to_core(model) = QIPresolve.build_model(QIPresolve.from_moi(backend(model)))

const RANDOM_QIP_DEFAULT_KWARGS = (
    p_con_eq = 0.5,
    var_threshold_lb = -3,
    var_threshold_ub = 4,
    p_var_is_candidate = 0.7,
    p_var_bilin = 0.4,
    p_var_diag = 0.5,
    p_var_lin = 0.6,
    coeff_lb = -5,
    coeff_ub = 5,
    force_diag_even = false,
    force_lin_even = false,
    force_feasibility = true,
    constraint_slack_range = -3:3,
)

random_qip_kwargs(; overrides...) = merge(RANDOM_QIP_DEFAULT_KWARGS, (; overrides...))

function expr_signature(expr::PC.QuadExpr)
    var_ids = sort!(collect(PC.vars(expr)))
    lin_terms = Tuple{Int, Float64}[]
    quad_terms = Tuple{Int, Int, Float64}[]

    for (idx, var_id_1) in enumerate(var_ids)
        lin_coeff = PC.get_lin_coeff(expr, var_id_1)
        isapprox(lin_coeff, 0.0) || push!(lin_terms, (var_id_1, lin_coeff))

        for var_id_2 in @view var_ids[idx:end]
            quad_coeff = PC.get_quad_coeff(expr, var_id_1, var_id_2)
            isapprox(quad_coeff, 0.0) || push!(quad_terms, (var_id_1, var_id_2, quad_coeff))
        end
    end

    return (
        constant = expr.constant,
        lin_terms = sort!(lin_terms),
        quad_terms = sort!(quad_terms),
    )
end

function qip_signature(model)
    qp_model = random_qip_to_core(model)
    return (
        vars = sort!([(var_id, var.lb, var.ub) for (var_id, var) in qp_model.vars]),
        cons = [
            (lhs = con.lhs, rhs = con.rhs, expr = expr_signature(con.qe))
            for con in qp_model.cons
        ],
        obj_sense = qp_model.obj_sense,
        obj_expr = expr_signature(qp_model.obj_expr),
    )
end

function nonzero_lin_terms(expr::PC.QuadExpr)
    return [
        (var_id, PC.get_lin_coeff(expr, var_id))
        for var_id in sort!(collect(PC.vars(expr)))
        if !isapprox(PC.get_lin_coeff(expr, var_id), 0.0)
    ]
end

function nonzero_quad_terms(expr::PC.QuadExpr)
    var_ids = sort!(collect(PC.vars(expr)))
    terms = Tuple{Int, Int, Float64}[]
    for (idx, var_id_1) in enumerate(var_ids)
        for var_id_2 in @view var_ids[idx:end]
            coeff = PC.get_quad_coeff(expr, var_id_1, var_id_2)
            isapprox(coeff, 0.0) || push!(terms, (var_id_1, var_id_2, coeff))
        end
    end
    return terms
end

is_even_coeff(coeff::Float64) = isinteger(coeff) && iseven(round(Int, coeff))

function remains_connected_after_removing_two(g::Graphs.AbstractGraph)
    n = nv(g)
    for first_vertex in 1:(n - 1)
        for second_vertex in (first_vertex + 1):n
            retained = [vertex for vertex in 1:n if vertex != first_vertex && vertex != second_vertex]
            subgraph, _ = induced_subgraph(g, retained)
            is_connected(subgraph) || return false
        end
    end
    return true
end

model_bounds(x, y) = [
    (
        JuMP.lower_bound(x[i]),
        JuMP.upper_bound(x[i]),
        JuMP.lower_bound(y[i]),
        JuMP.upper_bound(y[i]),
    )
    for i in eachindex(x)
]


@testset "graph instance public API" begin
    exported_names = names(IG)
    @test :generate_2_connected_instance in exported_names
    @test :generate_laman_instance in exported_names
    @test :generate_globally_rigid_instance in exported_names

    for helper in (
            :random_2_connected_graph,
            :plot_2_connected_graph,
            :random_laman_graph,
            :plot_laman_graph,
            :random_globally_rigid_graph,
            :plot_globally_rigid_graph,
        )
        @test isdefined(IG, helper)
        @test helper ∉ exported_names
    end
end


@testset "Laman graph generation" begin
    @test_throws ArgumentError random_laman_graph(2)
    @test_throws ArgumentError random_laman_graph(5; R = -1)
    @test_throws ArgumentError random_laman_graph(5; R = 0)
    @test_throws ArgumentError random_laman_graph(5; pH2 = -0.1)
    @test_throws ArgumentError random_laman_graph(5; pH2 = 1.1)
    @test_throws ArgumentError random_laman_graph(5; max_global_tries = 0)
    @test_throws ArgumentError random_laman_graph(5; max_tries_H1 = 0)
    @test_throws ArgumentError random_laman_graph(5; max_tries_H2 = 0)

    graph_1, coords_1 = random_laman_graph(8; R = 8, pH2 = 0.65, seed = 101)
    graph_2, coords_2 = random_laman_graph(8; R = 8, pH2 = 0.65, seed = 101)
    @test nv(graph_1) == 8
    @test ne(graph_1) == 2nv(graph_1) - 3
    @test canonical_edges(graph_1) == canonical_edges(graph_2)
    @test coord_tuples(coords_1) == coord_tuples(coords_2)
    @test allunique(coord_tuples(coords_1))
end


@testset "globally rigid graph generation" begin
    @test_throws ArgumentError random_globally_rigid_graph(3)
    @test_throws ArgumentError random_globally_rigid_graph(5; R = -1)
    @test_throws ArgumentError random_globally_rigid_graph(10; R = 1)
    @test_throws ArgumentError random_globally_rigid_graph(5; max_global_tries = 0)
    @test_throws ArgumentError random_globally_rigid_graph(5; max_tries_H2 = 0)

    base_graph, base_coords = random_globally_rigid_graph(4; R = 5, seed = 8)
    @test canonical_edges(base_graph) == canonical_edges(complete_graph(4))
    @test length(base_coords) == 4
    @test IG.no_three_collinear(base_coords)

    graph_1, coords_1 = random_globally_rigid_graph(9; R = 9, seed = 77)
    graph_2, coords_2 = random_globally_rigid_graph(9; R = 9, seed = 77)
    @test nv(graph_1) == 9
    @test ne(graph_1) == 2nv(graph_1) - 2
    @test canonical_edges(graph_1) == canonical_edges(graph_2)
    @test coord_tuples(coords_1) == coord_tuples(coords_2)
    @test allunique(coord_tuples(coords_1))
    @test remains_connected_after_removing_two(graph_1)
end


@testset "embedded model anchor bounds and symmetry" begin
    graph = complete_graph(4)
    coords = IG.IPoint[
        IG.IPoint(2, 3),
        IG.IPoint(5, 3),
        IG.IPoint(2, 7),
        IG.IPoint(6, 7),
    ]
    embedded = IG.to_embedded(graph, coords)

    center_refs, center_coords, center_distances = IG.embedding_references(embedded, Int[])
    @test center_refs == [2]
    @test center_coords == [IG.IPoint(0, 0)]
    @test center_distances[1, :] ≈ [3.0, 0.0, 5.0, sqrt(17.0)]

    model_0, x_0, y_0 = IG.build_embedding_model(embedded)
    @test model_bounds(x_0, y_0) == [
        (-3.0, 3.0, -3.0, 3.0),
        (0.0, 0.0, 0.0, 0.0),
        (-5.0, 5.0, -5.0, 5.0),
        (-4.0, 4.0, -4.0, 4.0),
    ]
    @test JuMP.num_constraints(model_0, JuMP.QuadExpr, MOI.EqualTo{Float64}) == 6
    @test JuMP.num_constraints(model_0, JuMP.QuadExpr, MOI.LessThan{Float64}) == 3
    @test JuMP.num_constraints(model_0, JuMP.AffExpr, MOI.EqualTo{Float64}) == 2
    @test JuMP.num_constraints(model_0, JuMP.AffExpr, MOI.LessThan{Float64}) == 1
    @test JuMP.num_constraints(model_0, JuMP.AffExpr, MOI.GreaterThan{Float64}) == 2

    anchor_refs, anchor_coords, anchor_distances = IG.embedding_references(embedded, [2])
    @test anchor_refs == [2]
    @test anchor_coords == [coords[2]]
    @test anchor_distances[1, :] ≈ [3.0, 0.0, 5.0, sqrt(17.0)]

    model_1, x_1, y_1 = IG.build_embedding_model(embedded, [2])
    @test model_bounds(x_1, y_1) == [
        (2.0, 8.0, 0.0, 6.0),
        (5.0, 5.0, 3.0, 3.0),
        (0.0, 10.0, -2.0, 8.0),
        (1.0, 9.0, -1.0, 7.0),
    ]
    @test JuMP.num_constraints(model_1, JuMP.QuadExpr, MOI.LessThan{Float64}) == 3
    @test JuMP.num_constraints(model_1, JuMP.AffExpr, MOI.EqualTo{Float64}) == 2
    @test JuMP.num_constraints(model_1, JuMP.AffExpr, MOI.LessThan{Float64}) == 1
    @test JuMP.num_constraints(model_1, JuMP.AffExpr, MOI.GreaterThan{Float64}) == 2
    diagonal_constraint = only(JuMP.all_constraints(model_1, JuMP.AffExpr, MOI.LessThan{Float64}))
    @test JuMP.constraint_object(diagonal_constraint).set.upper == -2.0
    quadrant_constraints = JuMP.all_constraints(model_1, JuMP.AffExpr, MOI.GreaterThan{Float64})
    @test sort([JuMP.constraint_object(constraint).set.lower for constraint in quadrant_constraints]) == [3.0, 5.0]

    reference_vertices, reference_coords, distances = IG.embedding_references(embedded, [1, 4])
    x_lower, x_upper, y_lower, y_upper = IG.merged_coordinate_bounds(reference_coords, distances)
    @test reference_vertices == [1, 4]
    @test reference_coords == coords[[1, 4]]
    @test x_lower == [2, 2, 2, 6]
    @test x_upper == [2, 5, 6, 6]
    @test y_lower == [3, 3, 3, 7]
    @test y_upper == [3, 6, 7, 7]

    model_2, x_2, y_2 = IG.build_embedding_model(embedded, [1, 4])
    @test model_bounds(x_2, y_2) == [
        (2.0, 2.0, 3.0, 3.0),
        (2.0, 5.0, 3.0, 6.0),
        (2.0, 6.0, 3.0, 7.0),
        (6.0, 6.0, 7.0, 7.0),
    ]
    @test JuMP.num_constraints(model_2, JuMP.QuadExpr, MOI.LessThan{Float64}) == 6
    @test JuMP.num_constraints(model_2, JuMP.AffExpr, MOI.EqualTo{Float64}) == 4
    @test JuMP.num_constraints(model_2, JuMP.AffExpr, MOI.LessThan{Float64}) == 0
    @test JuMP.num_constraints(model_2, JuMP.AffExpr, MOI.GreaterThan{Float64}) == 0

    @test_throws ArgumentError IG.build_embedding_model(embedded, [0])
    @test_throws ArgumentError IG.build_embedding_model(embedded, [1, 1])
    @test_throws ArgumentError IG.build_embedding_model(embedded, [5])
end


@testset "anchor counts across graph instance generators" begin
    generators = (
        generate_2_connected_instance,
        generate_laman_instance,
        generate_globally_rigid_instance,
    )

    for generator in generators
        for num_anchors in (0, 1, 2, 6)
            model, x, y = generator(6; R = 8, seed = 29, num_anchors = num_anchors)
            @test model !== nothing
            @test length(x) == 6
            @test length(y) == 6
        end
        @test_throws ArgumentError generator(6; num_anchors = -1)
        @test_throws ArgumentError generator(6; num_anchors = 7)

        _, x_1, y_1 = generator(6; R = 8, seed = 91, num_anchors = 2)
        _, x_2, y_2 = generator(6; R = 8, seed = 91, num_anchors = 2)
        @test model_bounds(x_1, y_1) == model_bounds(x_2, y_2)
    end
end


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

@testset "generate_random_qip_model validation" begin
    @test_throws ArgumentError generate_random_qip_model(
        0,
        1;
        random_qip_kwargs()...,
    )
    @test_throws ArgumentError generate_random_qip_model(
        2,
        -1;
        random_qip_kwargs()...,
    )
    @test_throws ArgumentError generate_random_qip_model(
        2,
        1;
        random_qip_kwargs(p_con_eq = 1.1)...,
    )
    @test_throws ArgumentError generate_random_qip_model(
        2,
        1;
        random_qip_kwargs(var_threshold_lb = 5, var_threshold_ub = 4)...,
    )
    @test_throws ArgumentError generate_random_qip_model(
        2,
        1;
        random_qip_kwargs(coeff_lb = 0, coeff_ub = 0)...,
    )
    @test_throws ArgumentError generate_random_qip_model(
        2,
        1;
        random_qip_kwargs(
            coeff_lb = -1,
            coeff_ub = 1,
            force_lin_even = true,
            p_var_lin = 1.0,
        )...,
    )
    @test_throws ArgumentError generate_random_qip_model(
        2,
        1;
        random_qip_kwargs(
            p_con_eq = 0.0,
            force_feasibility = true,
            constraint_slack_range = 1:3,
        )...,
    )
end

@testset "generate_random_qip_model is reproducible with seed" begin
    kwargs = random_qip_kwargs(
        p_con_eq = 0.4,
        p_var_is_candidate = 0.8,
        p_var_bilin = 0.5,
        p_var_diag = 0.5,
        p_var_lin = 0.5,
        seed = 44,
    )

    model_1, x_star_1 = generate_random_qip_model(5, 4; kwargs...)
    model_2, x_star_2 = generate_random_qip_model(5, 4; kwargs...)

    @test x_star_1 == x_star_2
    @test qip_signature(model_1) == qip_signature(model_2)
end

@testset "generate_random_qip_model forced feasibility" begin
    kwargs = random_qip_kwargs(
        p_con_eq = 0.0,
        p_var_is_candidate = 0.8,
        p_var_bilin = 0.5,
        p_var_diag = 0.5,
        p_var_lin = 0.5,
        force_feasibility = true,
        constraint_slack_range = -5:5,
        seed = 71,
    )

    model, x_star = generate_random_qip_model(6, 8; kwargs...)
    qp_model = random_qip_to_core(model)

    @test length(qp_model.vars) == 6
    @test length(qp_model.cons) == 8
    @test qp_model.obj_sense == :min

    for (var_id, var) in qp_model.vars
        @test kwargs.var_threshold_lb <= var.lb <= var.ub <= kwargs.var_threshold_ub
        @test var.lb <= x_star[var_id] <= var.ub
    end

    for con in qp_model.cons
        value = PC.eval_full(con.qe, x_star)
        @test con.lhs <= value <= con.rhs
    end
end

@testset "generate_random_qip_model can leave x_star infeasible" begin
    kwargs = random_qip_kwargs(
        p_con_eq = 1.0,
        p_var_is_candidate = 1.0,
        p_var_bilin = 0.0,
        p_var_diag = 1.0,
        p_var_lin = 0.0,
        force_feasibility = false,
        constraint_slack_range = 1:1,
        seed = 12,
    )

    model, x_star = generate_random_qip_model(4, 3; kwargs...)
    qp_model = random_qip_to_core(model)

    @test all(con -> con.lhs == con.rhs, qp_model.cons)
    @test all(con -> PC.eval_full(con.qe, x_star) < con.lhs, qp_model.cons)
end

@testset "generate_random_qip_model dense candidate coverage and objective" begin
    kwargs = random_qip_kwargs(
        p_con_eq = 1.0,
        var_threshold_lb = 0,
        var_threshold_ub = 3,
        p_var_is_candidate = 1.0,
        p_var_bilin = 1.0,
        p_var_diag = 1.0,
        p_var_lin = 1.0,
        coeff_lb = 1,
        coeff_ub = 1,
        force_feasibility = true,
        constraint_slack_range = 0:0,
        seed = 9,
    )

    model, _ = generate_random_qip_model(4, 1; kwargs...)
    qp_model = random_qip_to_core(model)
    con_expr = only(qp_model.cons).qe

    @test length(nonzero_quad_terms(con_expr)) == 10
    @test length(nonzero_lin_terms(con_expr)) == 4
    @test length(nonzero_lin_terms(qp_model.obj_expr)) == 4
    @test isempty(nonzero_quad_terms(qp_model.obj_expr))

    for var_id in 1:4
        @test PC.get_quad_coeff(con_expr, var_id, var_id) == 1.0
        @test PC.get_lin_coeff(con_expr, var_id) == 1.0
        @test PC.get_lin_coeff(qp_model.obj_expr, var_id) == 1.0
    end
    for first_id in 1:3
        for second_id in (first_id + 1):4
            @test PC.get_quad_coeff(con_expr, first_id, second_id) == 1.0
        end
    end
end

@testset "generate_random_qip_model respects forced even coefficients" begin
    kwargs = random_qip_kwargs(
        p_con_eq = 1.0,
        p_var_is_candidate = 1.0,
        p_var_bilin = 1.0,
        p_var_diag = 1.0,
        p_var_lin = 1.0,
        coeff_lb = -5,
        coeff_ub = 5,
        force_diag_even = true,
        force_lin_even = true,
        force_feasibility = true,
        constraint_slack_range = 0:0,
        seed = 29,
    )

    model, _ = generate_random_qip_model(4, 2; kwargs...)
    qp_model = random_qip_to_core(model)

    for con in qp_model.cons
        for (var_id_1, var_id_2, coeff) in nonzero_quad_terms(con.qe)
            var_id_1 == var_id_2 || continue
            @test is_even_coeff(coeff)
        end
        for (_, coeff) in nonzero_lin_terms(con.qe)
            @test is_even_coeff(coeff)
        end
    end

    for (_, coeff) in nonzero_lin_terms(qp_model.obj_expr)
        @test is_even_coeff(coeff)
    end
end
