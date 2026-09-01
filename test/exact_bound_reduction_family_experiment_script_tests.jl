using Test

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

@testset "exact bound reduction family experiment parses sweep lists" begin
    config = ExactBoundExperimentScript.build_config([
        "--count", "7",
        "--nvars", "3,5:6",
        "--domain-ubs", "1:3",
        "--max-distinct-coeffs", "1,3,10",
        "--domain-lb", "0",
    ])

    @test config.count == 7
    @test config.nvars == [3, 5, 6]
    @test config.domain_ubs == [1, 2, 3]
    @test config.max_distinct_coeffs == [1, 3, 10]
    @test config.domain_lb == 0

    alias_config = ExactBoundExperimentScript.build_config([
        "--domain-ub", "2,4",
    ])
    @test alias_config.domain_ubs == [2, 4]

    @test_throws ErrorException ExactBoundExperimentScript.build_config([
        "--max-distinct-coeffs", "all",
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

@testset "exact bound reduction family experiment sweeps cartesian product" begin
    config = ExactBoundExperimentScript.build_config([
        "--count", "1",
        "--nvars", "2,3",
        "--domain-ubs", "1:2",
        "--max-distinct-coeffs", "1,2",
        "--seed-base", "900",
    ])

    result = ExactBoundExperimentScript.run_experiment(config)
    expected = [
        (nvars, domain_ub, max_distinct_coeffs)
        for nvars in [2, 3]
        for domain_ub in [1, 2]
        for max_distinct_coeffs in [1, 2]
    ]

    @test length(result.rows) == 8
    @test [
        (row.nvars, row.domain_ub, row.max_distinct_coeffs)
        for row in result.rows
    ] == expected
    @test all(row -> row.constraints == 1, result.rows)
    @test all(row -> 0.0 <= row.opt_avg_red <= 1.0, result.rows)
    for key in expected
        @test result.generated_constraints[key] == 1
    end
end

@testset "exact bound reduction family experiment writes summaries" begin
    mktempdir() do dir
        output_path = joinpath(dir, "exact_bound_summary.csv")

        result = ExactBoundExperimentScript.main([
            "--count", "1",
            "--nvars", "2",
            "--domain-ubs", "1,2",
            "--max-distinct-coeffs", "1",
            "--seed-base", "777",
            "--output", output_path,
        ])

        @test result.config.count == 1
        @test result.config.nvars == [2]
        @test result.config.domain_ubs == [1, 2]
        @test result.config.max_distinct_coeffs == [1]
        @test length(result.rows) == 2
        @test isfile(output_path)

        header, rows = exact_bound_csv_table(output_path)
        @test header == [
            "nvars",
            "domain_lb",
            "domain_ub",
            "max_distinct_coeffs",
            "constraints",
            "exact_assignments_per_constraint",
            "total_optimal_relative_bound_reduction",
            "opt_avg_red",
            "avg_wall_time_sec_per_constraint",
        ]
        @test length(rows) == 2
        @test all(row -> length(row) == length(header), rows)
    end
end
