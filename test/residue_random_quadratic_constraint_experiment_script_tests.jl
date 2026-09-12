using Test

const RESIDUE_RANDOM_QUADRATIC_CONSTRAINT_EXPERIMENT_SCRIPT = joinpath(
    @__DIR__,
    "..",
    "scripts",
    "residue_random_quadratic_constraint_experiment.jl",
)

include(RESIDUE_RANDOM_QUADRATIC_CONSTRAINT_EXPERIMENT_SCRIPT)
const RandomQuadraticResidueScript = Main.ResidueRandomQuadraticConstraintExperiment

const EXPECTED_PRIMES_LT_128 = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97, 101, 103,
    107, 109, 113, 127,
]
const EXPECTED_PRIME_POWERS_LT_64 = [
    2, 3, 4, 5, 7, 8, 9, 11, 13,
    16, 17, 19, 23, 25, 27, 29, 31,
    32, 37, 41, 43, 47, 49, 53, 59, 61,
]
const EXPECTED_PRIME_POWERS_LT_128 = [
    2, 3, 4, 5, 7, 8, 9, 11, 13,
    16, 17, 19, 23, 25, 27, 29, 31,
    32, 37, 41, 43, 47, 49, 53, 59, 61,
    64, 67, 71, 73, 79, 81, 83, 89,
    97, 101, 103, 107, 109, 113, 121, 125, 127,
]

function random_quadratic_csv_table(path::AbstractString)
    lines = readlines(path)
    return split(first(lines), ','), split.(lines[2:end], ','; keepempty = true)
end

function random_quadratic_coefficients(con)
    terms = RandomQuadraticResidueScript.expression_terms(con.qe)
    coefficients = Int[]
    append!(coefficients, coefficient for (_, coefficient) in terms.lin_terms)
    append!(coefficients, coefficient for (_, _, coefficient) in terms.quad_terms)
    return coefficients
end

function random_quadratic_assert_same_sample(sample_a, sample_b)
    @test sample_a.x_star == sample_b.x_star
    @test sample_a.con.lhs == sample_b.con.lhs
    @test sample_a.con.rhs == sample_b.con.rhs
    @test sample_a.generation_tries == sample_b.generation_tries

    var_ids = sort!(collect(keys(sample_a.model.vars)))
    @test var_ids == sort!(collect(keys(sample_b.model.vars)))
    for (index, var_id) in enumerate(var_ids)
        @test sample_a.model.vars[var_id] == sample_b.model.vars[var_id]
        @test RandomQuadraticResidueScript.PC.get_lin_coeff(sample_a.con.qe, var_id) ==
            RandomQuadraticResidueScript.PC.get_lin_coeff(sample_b.con.qe, var_id)
        for other_id in @view var_ids[index:end]
            @test RandomQuadraticResidueScript.PC.get_quad_coeff(sample_a.con.qe, var_id, other_id) ==
                RandomQuadraticResidueScript.PC.get_quad_coeff(sample_b.con.qe, var_id, other_id)
        end
    end
end

@testset "random quadratic residue experiment uses primes below 128" begin
    strategies = RandomQuadraticResidueScript.strategy_specs()

    @test length(strategies) == 1
    @test only(strategies).name == "primes_lt_128"
    @test only(strategies).moduli == EXPECTED_PRIMES_LT_128
end

@testset "random quadratic residue experiment selects moduli strategies" begin
    prime_power_config = RandomQuadraticResidueScript.build_config([
        "--moduli-strategy", "prime_powers_lt_64",
    ])
    prime_power_strategy = only(RandomQuadraticResidueScript.strategy_specs(prime_power_config))

    @test prime_power_config.moduli_strategies == ["prime_powers_lt_64"]
    @test prime_power_strategy.name == "prime_powers_lt_64"
    @test prime_power_strategy.moduli == EXPECTED_PRIME_POWERS_LT_64

    prime_power_128_config = RandomQuadraticResidueScript.build_config([
        "--moduli-strategy", "prime_powers_lt_128",
    ])
    prime_power_128_strategy = only(
        RandomQuadraticResidueScript.strategy_specs(prime_power_128_config)
    )

    @test prime_power_128_config.moduli_strategies == ["prime_powers_lt_128"]
    @test prime_power_128_strategy.name == "prime_powers_lt_128"
    @test prime_power_128_strategy.moduli == EXPECTED_PRIME_POWERS_LT_128
    @test 64 in prime_power_128_strategy.moduli
    @test !(128 in prime_power_128_strategy.moduli)

    multi_config = RandomQuadraticResidueScript.build_config([
        "--count", "1",
        "--nvars", "3",
        "--domain-lb", "0",
        "--domain-ub", "2",
        "--density", "0.75",
        "--treewidth-threshold", "3",
        "--moduli-strategies", "primes_lt_128,prime_powers_lt_64,prime_powers_lt_128",
    ])
    result = RandomQuadraticResidueScript.run_experiment(multi_config)

    @test [row.strategy for row in result.rows] == [
        "primes_lt_128",
        "prime_powers_lt_64",
        "prime_powers_lt_128",
    ]
    @test [row.num_moduli for row in result.rows] == [
        length(EXPECTED_PRIMES_LT_128),
        length(EXPECTED_PRIME_POWERS_LT_64),
        length(EXPECTED_PRIME_POWERS_LT_128),
    ]
    @test all(row -> row.constraints == 1, result.rows)
    @test_throws ErrorException RandomQuadraticResidueScript.build_config([
        "--moduli-strategy", "not_a_strategy",
    ])
end

@testset "random quadratic residue experiment generation is deterministic" begin
    config = RandomQuadraticResidueScript.validate_config(
        RandomQuadraticResidueScript.CliConfig(
            count = 2,
            nvars = [4],
            seed_base = 123,
            domain_lb = 0,
            domain_ub = 2,
            density = 0.5,
            treewidth_threshold = 3,
        ),
    )

    sample_a = RandomQuadraticResidueScript.generate_constraint_sample(config, 4, 123)
    sample_b = RandomQuadraticResidueScript.generate_constraint_sample(config, 4, 123)
    random_quadratic_assert_same_sample(sample_a, sample_b)

    result = RandomQuadraticResidueScript.run_experiment(config)
    @test length(result.rows) == 1
    @test only(result.rows).constraints == 2
    @test result.generated_constraints[4] == 2
end

@testset "random quadratic residue experiment samples all coefficients at density one" begin
    config = RandomQuadraticResidueScript.build_config([
        "--count", "1",
        "--nvars", "3",
        "--domain-lb", "0",
        "--domain-ub", "2",
        "--density", "1.0",
        "--coeff-lb", "1",
        "--coeff-ub", "1",
        "--max-distinct-coeffs", "1",
    ])

    sample = RandomQuadraticResidueScript.generate_constraint_sample(config, 3, 456)
    terms = RandomQuadraticResidueScript.expression_terms(sample.con.qe)
    quad_pairs = sort!([(first_id, second_id) for (first_id, second_id, _) in terms.quad_terms])

    @test length(terms.lin_terms) == 3
    @test sort!([var_id for (var_id, _) in terms.lin_terms]) == [1, 2, 3]
    @test length(terms.quad_terms) == 6
    @test quad_pairs == [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]
end

@testset "random quadratic residue experiment respects max distinct coefficients" begin
    config = RandomQuadraticResidueScript.build_config([
        "--count", "1",
        "--nvars", "4",
        "--domain-lb", "0",
        "--domain-ub", "2",
        "--density", "1.0",
        "--coeff-lb", "-5",
        "--coeff-ub", "5",
        "--max-distinct-coeffs", "1",
    ])

    sample = RandomQuadraticResidueScript.generate_constraint_sample(config, 4, 789)
    coefficients = random_quadratic_coefficients(sample.con)

    @test !isempty(coefficients)
    @test length(unique(coefficients)) <= 1
end

@testset "random quadratic residue experiment computes range metrics" begin
    qe = RandomQuadraticResidueScript.PC.QuadExpr(
        RandomQuadraticResidueScript.QuadTerm[],
        RandomQuadraticResidueScript.LinTerm[],
    )
    before = (lhs = 1.0, rhs = 9.0)
    interval_con = RandomQuadraticResidueScript.PC.Constraint(1, qe, 3.0, 7.0)
    equality_con = RandomQuadraticResidueScript.PC.Constraint(2, qe, 5.0, 5.0)

    @test RandomQuadraticResidueScript.Common.relative_bound_range_reduction(
        before,
        interval_con,
    ) == 0.5
    @test RandomQuadraticResidueScript.Common.constraint_tightened_to_equality(
        before,
        equality_con,
    )
    @test !RandomQuadraticResidueScript.Common.constraint_tightened_to_equality(
        before,
        interval_con,
    )
end

@testset "random quadratic residue experiment writes CSV summaries" begin
    mktempdir() do dir
        output_path = joinpath(dir, "random_quadratic_residue_summary.csv")

        result = RandomQuadraticResidueScript.main([
            "--count", "2",
            "--nvars", "3",
            "--seed-base", "900",
            "--domain-lb", "0",
            "--domain-ub", "2",
            "--density", "0.75",
            "--treewidth-threshold", "3",
            "--exact-enumeration", "true",
            "--output", output_path,
        ])

        @test result.config.count == 2
        @test result.config.nvars == [3]
        @test result.config.exact_enumeration
        @test length(result.rows) == 1
        row = only(result.rows)
        @test row.strategy == "primes_lt_128"
        @test row.moduli == join(EXPECTED_PRIMES_LT_128, " ")
        @test row.num_moduli == length(EXPECTED_PRIMES_LT_128)
        @test row.constraints == 2
        @test row.bounds_considered == 4
        @test 0.0 <= row.pct_constraints_tightened_to_equality <= 100.0
        @test row.avg_relative_bound_range_reduction >= 0.0
        @test row.exact_assignments_per_constraint == 27
        @test 0.0 <= row.pct_bounds_fully_tightened_to_optimal <= 100.0
        @test 0.0 <= row.pct_bounds_tightened <= 100.0
        @test row.pct_bounds_fully_tightened_to_optimal <= row.pct_bounds_tightened
        @test row.avg_bound_gap_to_optimal >= 0.0
        @test row.total_residue_time_sec >= 0.0
        @test row.avg_wall_time_sec_per_constraint >= 0.0
        @test isfile(output_path)

        header, rows = random_quadratic_csv_table(output_path)
        @test header == [
            "nvars",
            "density",
            "domain_lb",
            "domain_ub",
            "max_distinct_coeffs",
            "exact_enumeration",
            "strategy",
            "moduli",
            "num_moduli",
            "constraints",
            "bounds_considered",
            "pct_constraints_tightened_to_equality",
            "avg_relative_bound_range_reduction",
            "exact_assignments_per_constraint",
            "bounds_fully_tightened_to_optimal",
            "pct_bounds_fully_tightened_to_optimal",
            "bounds_tightened",
            "pct_bounds_tightened",
            "avg_bound_gap_to_optimal",
            "total_residue_time_sec",
            "avg_wall_time_sec_per_constraint",
        ]
        @test length(rows) == 1
        @test all(row -> length(row) == length(header), rows)
    end
end

@testset "random quadratic residue experiment can skip exact enumeration" begin
    mktempdir() do dir
        output_path = joinpath(dir, "random_quadratic_residue_summary_no_exact.csv")

        result = RandomQuadraticResidueScript.main([
            "--count", "2",
            "--nvars", "3",
            "--seed-base", "901",
            "--domain-lb", "0",
            "--domain-ub", "2",
            "--density", "0.75",
            "--treewidth-threshold", "3",
            "--exact-enumeration", "false",
            "--output", output_path,
        ])

        row = only(result.rows)
        @test !result.config.exact_enumeration
        @test row.exact_enumeration == false
        @test row.constraints == 2
        @test row.bounds_considered == 4
        @test row.exact_assignments_per_constraint === missing
        @test row.bounds_fully_tightened_to_optimal === missing
        @test row.pct_bounds_fully_tightened_to_optimal === missing
        @test row.avg_bound_gap_to_optimal === missing
        @test 0.0 <= row.pct_constraints_tightened_to_equality <= 100.0
        @test row.avg_relative_bound_range_reduction >= 0.0
        @test row.bounds_tightened >= 0
        @test 0.0 <= row.pct_bounds_tightened <= 100.0
        @test isfile(output_path)
    end
end
