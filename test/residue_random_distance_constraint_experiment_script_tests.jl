using Test

const RESIDUE_RANDOM_DISTANCE_CONSTRAINT_EXPERIMENT_SCRIPT = joinpath(
    @__DIR__,
    "..",
    "scripts",
    "residue_random_distance_constraint_experiment.jl",
)

include(RESIDUE_RANDOM_DISTANCE_CONSTRAINT_EXPERIMENT_SCRIPT)
const RandomDistanceResidueScript = Main.ResidueRandomDistanceConstraintExperiment

const DISTANCE_EXPECTED_PRIMES_LT_128 = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97, 101, 103,
    107, 109, 113, 127,
]
const DISTANCE_EXPECTED_PRIMES_LT_64 = [
    prime for prime in DISTANCE_EXPECTED_PRIMES_LT_128 if prime < 64
]
const DISTANCE_EXPECTED_PRIME_POWERS_LT_64 = [
    2, 3, 4, 5, 7, 8, 9, 11, 13,
    16, 17, 19, 23, 25, 27, 29, 31,
    32, 37, 41, 43, 47, 49, 53, 59, 61,
]
const DISTANCE_EXPECTED_PRIME_POWERS_LT_128 = [
    2, 3, 4, 5, 7, 8, 9, 11, 13,
    16, 17, 19, 23, 25, 27, 29, 31,
    32, 37, 41, 43, 47, 49, 53, 59, 61,
    64, 67, 71, 73, 79, 81, 83, 89,
    97, 101, 103, 107, 109, 113, 121, 125, 127,
]
const DISTANCE_EXPECTED_ALL_2_TO_63 = collect(2:63)
const DISTANCE_EXPECTED_ALL_2_TO_128 = collect(2:128)
const DISTANCE_EXPECTED_MODULI_FAMILY_NAMES = [
    "primes_lt_64",
    "prime_powers_lt_64",
    "all_2_to_63",
    "primes_lt_128",
    "prime_powers_lt_128",
    "all_2_to_128",
]
const DISTANCE_EXPECTED_MODULI_FAMILY_VALUES = [
    DISTANCE_EXPECTED_PRIMES_LT_64,
    DISTANCE_EXPECTED_PRIME_POWERS_LT_64,
    DISTANCE_EXPECTED_ALL_2_TO_63,
    DISTANCE_EXPECTED_PRIMES_LT_128,
    DISTANCE_EXPECTED_PRIME_POWERS_LT_128,
    DISTANCE_EXPECTED_ALL_2_TO_128,
]

function random_distance_csv_table(path::AbstractString)
    lines = readlines(path)
    return split(first(lines), ','), split.(lines[2:end], ','; keepempty = true)
end

function random_distance_assert_same_sample(sample_a, sample_b)
    @test sample_a.point_1 == sample_b.point_1
    @test sample_a.point_2 == sample_b.point_2
    @test sample_a.squared_distance == sample_b.squared_distance
    @test sample_a.x_star == sample_b.x_star
    @test sample_a.con.lhs == sample_b.con.lhs
    @test sample_a.con.rhs == sample_b.con.rhs
    @test sample_a.generation_tries == sample_b.generation_tries

    var_ids = sort!(collect(keys(sample_a.model.vars)))
    @test var_ids == [1, 2, 3, 4]
    @test var_ids == sort!(collect(keys(sample_b.model.vars)))
    for (index, var_id) in enumerate(var_ids)
        @test sample_a.model.vars[var_id] == sample_b.model.vars[var_id]
        @test RandomDistanceResidueScript.PC.get_lin_coeff(sample_a.con.qe, var_id) ==
            RandomDistanceResidueScript.PC.get_lin_coeff(sample_b.con.qe, var_id)
        for other_id in @view var_ids[index:end]
            @test RandomDistanceResidueScript.PC.get_quad_coeff(sample_a.con.qe, var_id, other_id) ==
                RandomDistanceResidueScript.PC.get_quad_coeff(sample_b.con.qe, var_id, other_id)
        end
    end
end

@testset "random distance residue experiment uses moduli family by default" begin
    strategies = RandomDistanceResidueScript.strategy_specs()

    @test length(strategies) == 6
    @test [strategy.name for strategy in strategies] == DISTANCE_EXPECTED_MODULI_FAMILY_NAMES
    @test [strategy.moduli for strategy in strategies] == DISTANCE_EXPECTED_MODULI_FAMILY_VALUES
end

@testset "random distance residue experiment selects moduli strategies" begin
    prime_power_config = RandomDistanceResidueScript.build_config([
        "--moduli-strategy", "prime_powers_lt_64",
    ])
    prime_power_strategy = only(RandomDistanceResidueScript.strategy_specs(prime_power_config))

    @test prime_power_config.moduli_strategies == ["prime_powers_lt_64"]
    @test prime_power_strategy.name == "prime_powers_lt_64"
    @test prime_power_strategy.moduli == DISTANCE_EXPECTED_PRIME_POWERS_LT_64

    prime_power_128_config = RandomDistanceResidueScript.build_config([
        "--moduli-strategy", "prime_powers_lt_128",
    ])
    prime_power_128_strategy = only(
        RandomDistanceResidueScript.strategy_specs(prime_power_128_config)
    )

    @test prime_power_128_config.moduli_strategies == ["prime_powers_lt_128"]
    @test prime_power_128_strategy.name == "prime_powers_lt_128"
    @test prime_power_128_strategy.moduli == DISTANCE_EXPECTED_PRIME_POWERS_LT_128
    @test 64 in prime_power_128_strategy.moduli
    @test !(128 in prime_power_128_strategy.moduli)

    all_63_config = RandomDistanceResidueScript.build_config([
        "--moduli-strategy", "all_2_to_63",
    ])
    all_63_strategy = only(RandomDistanceResidueScript.strategy_specs(all_63_config))

    @test all_63_config.moduli_strategies == ["all_2_to_63"]
    @test all_63_strategy.name == "all_2_to_63"
    @test all_63_strategy.moduli == DISTANCE_EXPECTED_ALL_2_TO_63

    all_128_config = RandomDistanceResidueScript.build_config([
        "--moduli-strategy", "all_2_to_128",
    ])
    all_128_strategy = only(RandomDistanceResidueScript.strategy_specs(all_128_config))

    @test all_128_config.moduli_strategies == ["all_2_to_128"]
    @test all_128_strategy.name == "all_2_to_128"
    @test all_128_strategy.moduli == DISTANCE_EXPECTED_ALL_2_TO_128

    multi_config = RandomDistanceResidueScript.build_config([
        "--count", "1",
        "--R", "2",
        "--alpha", "0.2",
        "--treewidth-threshold", "3",
        "--moduli-strategies", "primes_lt_128,prime_powers_lt_64,prime_powers_lt_128",
    ])
    result = RandomDistanceResidueScript.run_experiment(multi_config)

    @test [row.strategy for row in result.rows] == [
        "primes_lt_128",
        "prime_powers_lt_64",
        "prime_powers_lt_128",
    ]
    @test [row.num_moduli for row in result.rows] == [
        length(DISTANCE_EXPECTED_PRIMES_LT_128),
        length(DISTANCE_EXPECTED_PRIME_POWERS_LT_64),
        length(DISTANCE_EXPECTED_PRIME_POWERS_LT_128),
    ]
    @test all(row -> row.constraints == 1, result.rows)
    @test_throws ErrorException RandomDistanceResidueScript.build_config([
        "--moduli-strategy", "not_a_strategy",
    ])
end

@testset "random distance residue experiment generation is deterministic" begin
    config = RandomDistanceResidueScript.validate_config(
        RandomDistanceResidueScript.CliConfig(
            count = 2,
            R = 2,
            alpha = 0.2,
            seed_base = 123,
            treewidth_threshold = 3,
        ),
    )

    sample_a = RandomDistanceResidueScript.generate_constraint_sample(config, 123)
    sample_b = RandomDistanceResidueScript.generate_constraint_sample(config, 123)
    random_distance_assert_same_sample(sample_a, sample_b)

    result = RandomDistanceResidueScript.run_experiment(config)
    @test length(result.rows) == length(DISTANCE_EXPECTED_MODULI_FAMILY_NAMES)
    @test [row.strategy for row in result.rows] == DISTANCE_EXPECTED_MODULI_FAMILY_NAMES
    @test all(row -> row.constraints == 2, result.rows)
    @test result.generated_constraints == 2
end

@testset "random distance residue experiment builds squared distance intervals" begin
    config = RandomDistanceResidueScript.build_config([
        "--count", "1",
        "--R", "3",
        "--alpha", "0.25",
        "--seed-base", "456",
    ])

    sample = RandomDistanceResidueScript.generate_constraint_sample(config, 456)
    point_1 = sample.point_1
    point_2 = sample.point_2

    @test point_1 != point_2
    @test -config.R <= point_1.x <= config.R
    @test -config.R <= point_1.y <= config.R
    @test -config.R <= point_2.x <= config.R
    @test -config.R <= point_2.y <= config.R

    d2 = (point_1.x - point_2.x)^2 + (point_1.y - point_2.y)^2
    @test sample.squared_distance == d2
    @test RandomDistanceResidueScript.PC.eval_full(sample.con.qe, sample.x_star) == d2
    @test sample.con.lhs == ceil((1.0 - config.alpha) * d2)
    @test sample.con.rhs == floor((1.0 + config.alpha) * d2)

    @test RandomDistanceResidueScript.PC.get_quad_coeff(sample.con.qe, 1, 1) == 1.0
    @test RandomDistanceResidueScript.PC.get_quad_coeff(sample.con.qe, 2, 2) == 1.0
    @test RandomDistanceResidueScript.PC.get_quad_coeff(sample.con.qe, 3, 3) == 1.0
    @test RandomDistanceResidueScript.PC.get_quad_coeff(sample.con.qe, 4, 4) == 1.0
    @test RandomDistanceResidueScript.PC.get_quad_coeff(sample.con.qe, 1, 3) == -2.0
    @test RandomDistanceResidueScript.PC.get_quad_coeff(sample.con.qe, 2, 4) == -2.0
end

@testset "random distance residue experiment exact enumeration uses full grid domains" begin
    config = RandomDistanceResidueScript.build_config([
        "--count", "1",
        "--R", "2",
        "--alpha", "0.2",
        "--seed-base", "789",
    ])

    sample = RandomDistanceResidueScript.generate_constraint_sample(config, 789)
    exact = RandomDistanceResidueScript.exact_bound_tightening(sample.con, sample.model.vars)

    @test exact.assignment_count == (2 * config.R + 1)^4
    @test exact.lhs >= sample.con.lhs
    @test exact.rhs <= sample.con.rhs
end

@testset "random distance residue experiment writes CSV summaries" begin
    mktempdir() do dir
        output_path = joinpath(dir, "random_distance_residue_summary.csv")

        result = RandomDistanceResidueScript.main([
            "--count", "2",
            "--R", "2",
            "--alpha", "0.2",
            "--seed-base", "900",
            "--treewidth-threshold", "3",
            "--exact-enumeration", "true",
            "--output", output_path,
        ])

        @test result.config.count == 2
        @test result.config.R == 2
        @test result.config.alpha == 0.2
        @test result.config.exact_enumeration
        @test length(result.rows) == length(DISTANCE_EXPECTED_MODULI_FAMILY_NAMES)
        @test [row.strategy for row in result.rows] == DISTANCE_EXPECTED_MODULI_FAMILY_NAMES
        for (row, expected_moduli) in zip(result.rows, DISTANCE_EXPECTED_MODULI_FAMILY_VALUES)
            @test row.R == 2
            @test row.alpha == 0.2
            @test row.nvars == 4
            @test row.domain_lb == -2
            @test row.domain_ub == 2
            @test row.moduli == join(expected_moduli, " ")
            @test row.num_moduli == length(expected_moduli)
            @test row.constraints == 2
            @test row.bounds_considered == 4
            @test 0.0 <= row.pct_constraints_tightened_to_equality <= 100.0
            @test row.avg_relative_bound_range_reduction >= 0.0
            @test row.exact_assignments_per_constraint == 625
            @test 0.0 <= row.pct_bounds_fully_tightened_to_optimal <= 100.0
            @test 0.0 <= row.pct_bounds_tightened <= 100.0
            @test row.pct_bounds_fully_tightened_to_optimal <= row.pct_bounds_tightened
            @test row.avg_bound_gap_to_optimal >= 0.0
            @test row.total_residue_time_sec >= 0.0
            @test row.avg_wall_time_sec_per_constraint >= 0.0
        end
        @test isfile(output_path)

        header, rows = random_distance_csv_table(output_path)
        @test header == [
            "R",
            "alpha",
            "nvars",
            "domain_lb",
            "domain_ub",
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
        @test length(rows) == length(DISTANCE_EXPECTED_MODULI_FAMILY_NAMES)
        @test all(row -> length(row) == length(header), rows)
    end
end

@testset "random distance residue experiment can skip exact enumeration" begin
    mktempdir() do dir
        output_path = joinpath(dir, "random_distance_residue_summary_no_exact.csv")

        result = RandomDistanceResidueScript.main([
            "--count", "2",
            "--R", "2",
            "--alpha", "0.2",
            "--seed-base", "901",
            "--treewidth-threshold", "3",
            "--exact-enumeration", "false",
            "--output", output_path,
        ])

        @test !result.config.exact_enumeration
        @test length(result.rows) == length(DISTANCE_EXPECTED_MODULI_FAMILY_NAMES)
        @test [row.strategy for row in result.rows] == DISTANCE_EXPECTED_MODULI_FAMILY_NAMES
        for row in result.rows
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
        end
        @test isfile(output_path)
    end
end
