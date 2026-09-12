using Test
using Statistics

const EXACT_BOUND_REDUCTION_FAMILY_EXPERIMENT_SCRIPT = joinpath(
    @__DIR__,
    "..",
    "scripts",
    "exact_bound_reduction_family_experiment.jl",
)

include(EXACT_BOUND_REDUCTION_FAMILY_EXPERIMENT_SCRIPT)
const ExactBoundExperimentScript = Main.ExactBoundReductionFamilyExperiment

function exact_bound_csv_table(path::AbstractString)
    lines = readlines(path)
    return split(first(lines), ','), split.(lines[2:end], ','; keepempty = true)
end

function exact_bound_treewidth_constraint(nvars::Int, edges)
    quad_terms = Tuple{Float64, ExactBoundExperimentScript.PC.VarId, ExactBoundExperimentScript.PC.VarId}[
        (2.0, first_var, second_var)
        for (first_var, second_var) in edges
    ]
    lin_terms = Tuple{Float64, ExactBoundExperimentScript.PC.VarId}[]
    return ExactBoundExperimentScript.PC.Constraint(
        1,
        ExactBoundExperimentScript.PC.QuadExpr(quad_terms, lin_terms),
        0.0,
        1.0,
    )
end

function complete_graph_edges(nvars::Int)
    return [
        (first_var, second_var)
        for first_var in 1:(nvars - 1)
        for second_var in (first_var + 1):nvars
    ]
end

@testset "exact bound reduction family experiment parses sweep lists" begin
    config = ExactBoundExperimentScript.build_config([
        "--count", "7",
        "--nvars", "3,5:6",
        "--density", "0.1,0.35",
        "--domain-ubs", "1:3",
        "--max-distinct-coeffs", "1,3,10",
        "--domain-lb", "0",
    ])

    @test config.count == 7
    @test config.nvars == [3, 5, 6]
    @test config.densities == [0.1, 0.35]
    @test config.domain_ubs == [1, 2, 3]
    @test config.max_distinct_coeffs == [1, 3, 10]
    @test config.domain_lb == 0

    alias_config = ExactBoundExperimentScript.build_config([
        "--domain-ub", "2,4",
    ])
    @test alias_config.domain_ubs == [2, 4]

    output_config = ExactBoundExperimentScript.build_config([
        "--output-dir", "exact_output",
    ])
    @test basename(output_config.output_dir) == "exact_output"

    @test_throws ErrorException ExactBoundExperimentScript.build_config([
        "--max-distinct-coeffs", "all",
    ])
    @test_throws ErrorException ExactBoundExperimentScript.build_config([
        "--density", "0.1,1.2",
    ])
    @test_throws ErrorException ExactBoundExperimentScript.build_config([
        "--diag-probability", "0.1",
    ])
end

@testset "exact bound reduction family experiment generation is deterministic" begin
    config = ExactBoundExperimentScript.validate_config(
        ExactBoundExperimentScript.CliConfig(
            count = 2,
            nvars = [5],
            domain_ubs = [2],
            max_distinct_coeffs = [3],
            seed_base = 123,
        ),
    )

    sample_a = ExactBoundExperimentScript.generate_constraint_sample(config, 5, 2, 3, 123)
    sample_b = ExactBoundExperimentScript.generate_constraint_sample(config, 5, 2, 3, 123)

    @test sample_a.x_star == sample_b.x_star
    @test sample_a.con.lhs == sample_b.con.lhs
    @test sample_a.con.rhs == sample_b.con.rhs
    @test sample_a.domain_ub == sample_b.domain_ub == 2
    @test sample_a.max_distinct_coeffs == sample_b.max_distinct_coeffs == 3

    var_ids = sort!(collect(keys(sample_a.model.vars)))
    for (index, var_id) in enumerate(var_ids)
        @test sample_a.model.vars[var_id] == sample_b.model.vars[var_id]
        @test ExactBoundExperimentScript.PC.get_lin_coeff(sample_a.con.qe, var_id) ==
            ExactBoundExperimentScript.PC.get_lin_coeff(sample_b.con.qe, var_id)
        for other_id in @view var_ids[index:end]
            @test ExactBoundExperimentScript.PC.get_quad_coeff(sample_a.con.qe, var_id, other_id) ==
                ExactBoundExperimentScript.PC.get_quad_coeff(sample_b.con.qe, var_id, other_id)
        end
    end
end

@testset "exact bound reduction family experiment exact enumeration" begin
    vars = Dict{ExactBoundExperimentScript.PC.VarId, ExactBoundExperimentScript.PC.IntVar}(
        1 => ExactBoundExperimentScript.PC.IntVar(0.0, 2.0),
        2 => ExactBoundExperimentScript.PC.IntVar(0.0, 2.0),
    )
    con = ExactBoundExperimentScript.PC.Constraint(
        1,
        ExactBoundExperimentScript.PC.QuadExpr(
            Tuple{Float64, ExactBoundExperimentScript.PC.VarId, ExactBoundExperimentScript.PC.VarId}[
                (2.0, 1, 2),
            ],
            Tuple{Float64, ExactBoundExperimentScript.PC.VarId}[
                (1.0, 1),
            ],
        ),
        4.0,
        9.0,
    )

    exact = ExactBoundExperimentScript.exact_bound_tightening(con, vars)

    @test exact.lhs == 5.0
    @test exact.rhs == 6.0
    @test exact.assignment_count == 9
    @test exact.relative_bound_reduction == 0.8
end

@testset "exact bound reduction family experiment summarizes reduction variance" begin
    result = ExactBoundExperimentScript.SweepResult(
        nvars = 2,
        density = 0.1,
        domain_lb = 0,
        domain_ub = 1,
        max_distinct_coeffs = 1,
        constraints = 3,
        optimal_relative_bound_reductions = [0.0, 0.5, 1.0],
        treewidths = [1, 2, 4],
    )

    row = ExactBoundExperimentScript.result_row(result)

    @test row.opt_avg_red == mean(result.optimal_relative_bound_reductions)
    @test row.opt_std_red == std(result.optimal_relative_bound_reductions)
    @test row.tw_1 ≈ 100 / 3
    @test row.tw_2 ≈ 100 / 3
    @test row.tw_3 == 0.0
    @test row.tw_ge4 ≈ 100 / 3

    single_result = ExactBoundExperimentScript.SweepResult(
        nvars = 2,
        density = 0.1,
        domain_lb = 0,
        domain_ub = 1,
        max_distinct_coeffs = 1,
        constraints = 1,
        optimal_relative_bound_reductions = [0.5],
        treewidths = [0],
    )

    @test isnan(ExactBoundExperimentScript.result_row(single_result).opt_std_red)
end

@testset "exact bound reduction family experiment computes interaction treewidth" begin
    path = exact_bound_treewidth_constraint(4, [(1, 2), (2, 3), (3, 4)])
    triangle = exact_bound_treewidth_constraint(3, complete_graph_edges(3))
    complete_four = exact_bound_treewidth_constraint(4, complete_graph_edges(4))
    complete_five = exact_bound_treewidth_constraint(5, complete_graph_edges(5))
    singleton = ExactBoundExperimentScript.PC.Constraint(
        1,
        ExactBoundExperimentScript.PC.QuadExpr(
            Tuple{Float64, ExactBoundExperimentScript.PC.VarId, ExactBoundExperimentScript.PC.VarId}[],
            Tuple{Float64, ExactBoundExperimentScript.PC.VarId}[(1.0, 1)],
        ),
        0.0,
        1.0,
    )

    @test ExactBoundExperimentScript.constraint_treewidth(path) == 1
    @test ExactBoundExperimentScript.constraint_treewidth(triangle) == 2
    @test ExactBoundExperimentScript.constraint_treewidth(complete_four) == 3
    @test ExactBoundExperimentScript.constraint_treewidth(complete_five) >= 4
    @test ExactBoundExperimentScript.constraint_treewidth(singleton) == 0
end

@testset "exact bound reduction family experiment aggregates grouped rows" begin
    results = [
        ExactBoundExperimentScript.SweepResult(
            nvars = 2,
            density = 0.1,
            domain_lb = 0,
            domain_ub = 1,
            max_distinct_coeffs = 3,
            constraints = 2,
            optimal_relative_bound_reductions = [0.0, 0.5],
            treewidths = [1, 2],
            total_exact_time_sec = 2.0,
        ),
        ExactBoundExperimentScript.SweepResult(
            nvars = 2,
            density = 0.2,
            domain_lb = 0,
            domain_ub = 1,
            max_distinct_coeffs = 5,
            constraints = 1,
            optimal_relative_bound_reductions = [1.0],
            treewidths = [4],
            total_exact_time_sec = 4.0,
        ),
        ExactBoundExperimentScript.SweepResult(
            nvars = 3,
            density = 0.1,
            domain_lb = 0,
            domain_ub = 2,
            max_distinct_coeffs = 3,
            constraints = 2,
            optimal_relative_bound_reductions = [0.25, 0.75],
            treewidths = [3, 5],
            total_exact_time_sec = 6.0,
        ),
    ]

    aggregate_rows = ExactBoundExperimentScript.aggregate_result_rows(results)
    nvars_two = only(row for row in aggregate_rows.by_nvars if row.nvars == 2)
    density_point_one = only(row for row in aggregate_rows.by_density if row.density == 0.1)
    domain_one = only(row for row in aggregate_rows.by_domain_ub if row.domain_ub == 1)
    coeffs_three = only(
        row for row in aggregate_rows.by_max_distinct_coeffs if row.max_distinct_coeffs == 3
    )

    @test nvars_two.constraints == 3
    @test nvars_two.opt_avg_red == mean([0.0, 0.5, 1.0])
    @test nvars_two.opt_std_red == std([0.0, 0.5, 1.0])
    @test nvars_two.avg_wall_time_sec_per_constraint == 2.0
    @test nvars_two.tw_1 ≈ 100 / 3
    @test nvars_two.tw_2 ≈ 100 / 3
    @test nvars_two.tw_3 == 0.0
    @test nvars_two.tw_ge4 ≈ 100 / 3

    @test density_point_one.constraints == 4
    @test density_point_one.opt_avg_red == mean([0.0, 0.5, 0.25, 0.75])
    @test density_point_one.opt_std_red == std([0.0, 0.5, 0.25, 0.75])
    @test density_point_one.avg_wall_time_sec_per_constraint == 2.0
    @test density_point_one.tw_1 == 25.0
    @test density_point_one.tw_2 == 25.0
    @test density_point_one.tw_3 == 25.0
    @test density_point_one.tw_ge4 == 25.0

    @test domain_one.constraints == 3
    @test coeffs_three.constraints == 4
end

@testset "exact bound reduction family experiment sweeps cartesian product" begin
    config = ExactBoundExperimentScript.build_config([
        "--count", "1",
        "--nvars", "2,3",
        "--densities", "0.1,0.2",
        "--domain-ubs", "1:2",
        "--max-distinct-coeffs", "1,2",
        "--seed-base", "900",
    ])

    result = ExactBoundExperimentScript.run_experiment(config)
    expected = [
        (nvars, density, domain_ub, max_distinct_coeffs)
        for nvars in [2, 3]
        for density in [0.1, 0.2]
        for domain_ub in [1, 2]
        for max_distinct_coeffs in [1, 2]
    ]

    @test length(result.rows) == 16
    @test [
        (row.nvars, row.density, row.domain_ub, row.max_distinct_coeffs)
        for row in result.rows
    ] == expected
    @test all(row -> row.constraints == 1, result.rows)
    @test all(row -> 0.0 <= row.opt_avg_red <= 1.0, result.rows)
    @test all(row -> isnan(row.opt_std_red), result.rows)
    @test all(row -> 0.0 <= row.tw_1 <= 100.0, result.rows)
    @test all(row -> 0.0 <= row.tw_2 <= 100.0, result.rows)
    @test all(row -> 0.0 <= row.tw_3 <= 100.0, result.rows)
    @test all(row -> 0.0 <= row.tw_ge4 <= 100.0, result.rows)
    @test all(row -> row.constraints == 8, result.aggregate_rows.by_nvars)
    @test all(row -> row.constraints == 8, result.aggregate_rows.by_density)
    @test all(row -> row.constraints == 8, result.aggregate_rows.by_domain_ub)
    @test all(row -> row.constraints == 8, result.aggregate_rows.by_max_distinct_coeffs)
    for key in expected
        @test result.generated_constraints[key] == 1
    end
end

@testset "exact bound reduction family experiment writes summaries" begin
    mktempdir() do dir
        output_dir = joinpath(dir, "exact_output")

        result = ExactBoundExperimentScript.main([
            "--count", "1",
            "--nvars", "2",
            "--domain-ubs", "1,2",
            "--max-distinct-coeffs", "1",
            "--seed-base", "777",
            "--output", output_dir,
        ])

        @test result.config.count == 1
        @test result.config.nvars == [2]
        @test result.config.densities == [0.1]
        @test result.config.domain_ubs == [1, 2]
        @test result.config.max_distinct_coeffs == [1]
        @test length(result.rows) == 2
        @test isdir(output_dir)

        expected_outputs = [
            (
                path = result.output_paths.by_nvars,
                field = "nvars",
                row_count = 1,
            ),
            (
                path = result.output_paths.by_density,
                field = "density",
                row_count = 1,
            ),
            (
                path = result.output_paths.by_domain_ub,
                field = "domain_ub",
                row_count = 2,
            ),
            (
                path = result.output_paths.by_max_distinct_coeffs,
                field = "max_distinct_coeffs",
                row_count = 1,
            ),
        ]

        for expected in expected_outputs
            @test isfile(expected.path)
            header, rows = exact_bound_csv_table(expected.path)
            @test header == [
                expected.field,
                "constraints",
                "opt_avg_red",
                "opt_std_red",
                "tw_1",
                "tw_2",
                "tw_3",
                "tw_ge4",
                "avg_wall_time_sec_per_constraint",
            ]
            @test length(rows) == expected.row_count
            @test all(row -> length(row) == length(header), rows)
        end
    end
end
