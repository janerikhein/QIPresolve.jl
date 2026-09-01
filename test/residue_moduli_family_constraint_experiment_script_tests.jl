using Test

const RESIDUE_MODULI_FAMILY_CONSTRAINT_EXPERIMENT_SCRIPT = joinpath(
    @__DIR__,
    "..",
    "scripts",
    "residue_moduli_family_constraint_experiment.jl",
)

include(RESIDUE_MODULI_FAMILY_CONSTRAINT_EXPERIMENT_SCRIPT)
const ResidueFamilyExperimentScript = Main.ResidueModuliFamilyConstraintExperiment

function residue_family_csv_table(path::AbstractString)
    lines = readlines(path)
    return split(first(lines), ','), split.(lines[2:end], ','; keepempty = true)
end

function residue_family_strategy_map(limit::Int = 64)
    return Dict(spec.name => spec.moduli for spec in ResidueFamilyExperimentScript.strategy_specs(limit))
end

function nonzero_coefficients(con)
    terms = ResidueFamilyExperimentScript.expression_terms(con.qe)
    coefficients = Int[]
    append!(coefficients, coefficient for (_, coefficient) in terms.lin_terms)
    append!(coefficients, coefficient for (_, _, coefficient) in terms.quad_terms)
    return coefficients
end

@testset "residue moduli family experiment builds requested families" begin
    strategies = residue_family_strategy_map()

    @test strategies["primes"] == [
        2, 3, 5, 7, 11, 13, 17, 19, 23,
        29, 31, 37, 41, 43, 47, 53, 59, 61,
    ]
    @test strategies["prime_powers"] == [
        2, 3, 4, 5, 7, 8, 9, 11, 13,
        16, 17, 19, 23, 25, 27, 29, 31,
        32, 37, 41, 43, 47, 49, 53, 59, 61,
    ]
    @test !(64 in strategies["prime_powers"])
    @test strategies["full"] == collect(2:63)
end

@testset "residue moduli family experiment generation is deterministic" begin
    config = ResidueFamilyExperimentScript.validate_config(
        ResidueFamilyExperimentScript.CliConfig(
            count = 2,
            nvars = [5],
            seed_base = 123,
            domain_lb = 0,
            domain_ub = 3,
            treewidth_threshold = 3,
        ),
    )

    sample_a = ResidueFamilyExperimentScript.generate_constraint_sample(config, 5, 123)
    sample_b = ResidueFamilyExperimentScript.generate_constraint_sample(config, 5, 123)

    @test sample_a.x_star == sample_b.x_star
    @test sample_a.con.lhs == sample_b.con.lhs
    @test sample_a.con.rhs == sample_b.con.rhs

    var_ids = sort!(collect(keys(sample_a.model.vars)))
    for (index, var_id) in enumerate(var_ids)
        @test sample_a.model.vars[var_id] == sample_b.model.vars[var_id]
        @test ResidueFamilyExperimentScript.PC.get_lin_coeff(sample_a.con.qe, var_id) ==
            ResidueFamilyExperimentScript.PC.get_lin_coeff(sample_b.con.qe, var_id)
        for other_id in @view var_ids[index:end]
            @test ResidueFamilyExperimentScript.PC.get_quad_coeff(sample_a.con.qe, var_id, other_id) ==
                ResidueFamilyExperimentScript.PC.get_quad_coeff(sample_b.con.qe, var_id, other_id)
        end
    end

    result = ResidueFamilyExperimentScript.run_experiment(config)
    @test length(result.rows) == 3
    @test all(row -> row.constraints == 2, result.rows)
    @test result.generated_constraints[5] == 2
end

@testset "residue moduli family experiment respects max distinct coefficients" begin
    config = ResidueFamilyExperimentScript.build_config([
        "--count", "1",
        "--nvars", "5",
        "--domain-lb", "0",
        "--domain-ub", "2",
        "--diag-probability", "1.0",
        "--linear-probability", "1.0",
        "--max-distinct-coeffs", "1",
    ])

    @test config.max_distinct_coeffs == 1

    sample = ResidueFamilyExperimentScript.generate_constraint_sample(config, 5, 456)
    coefficients = nonzero_coefficients(sample.con)

    @test !isempty(coefficients)
    @test length(unique(coefficients)) <= 1
end

@testset "residue moduli family experiment gcd scales generated constraints" begin
    config = ResidueFamilyExperimentScript.build_config([
        "--count", "1",
        "--nvars", "2",
        "--domain-lb", "0",
        "--domain-ub", "2",
        "--extra-edge-probability", "0.0",
        "--coeff-lb", "4",
        "--coeff-ub", "4",
        "--diag-probability", "0.0",
        "--linear-probability", "0.0",
    ])

    sample = ResidueFamilyExperimentScript.generate_constraint_sample(config, 2, 789)

    @test ResidueFamilyExperimentScript.PC.get_quad_coeff(sample.con.qe, 1, 2) == 2.0
end

@testset "residue moduli family experiment exact enumeration" begin
    vars = Dict{ResidueFamilyExperimentScript.PC.VarId, ResidueFamilyExperimentScript.PC.IntVar}(
        1 => ResidueFamilyExperimentScript.PC.IntVar(0.0, 2.0),
        2 => ResidueFamilyExperimentScript.PC.IntVar(0.0, 2.0),
    )
    con = ResidueFamilyExperimentScript.PC.Constraint(
        1,
        ResidueFamilyExperimentScript.PC.QuadExpr(
            Tuple{Float64, ResidueFamilyExperimentScript.PC.VarId, ResidueFamilyExperimentScript.PC.VarId}[
                (2.0, 1, 2),
            ],
            Tuple{Float64, ResidueFamilyExperimentScript.PC.VarId}[
                (1.0, 1),
            ],
        ),
        4.0,
        9.0,
    )

    exact = ResidueFamilyExperimentScript.exact_bound_tightening(con, vars)

    @test exact.lhs == 5.0
    @test exact.rhs == 6.0
    @test exact.assignment_count == 9
    @test exact.relative_bound_reduction == 0.8
end

@testset "residue moduli family experiment writes base summaries" begin
    mktempdir() do dir
        output_path = joinpath(dir, "residue_family_summary.csv")

        result = ResidueFamilyExperimentScript.main([
            "--count", "2",
            "--nvars", "5",
            "--seed-base", "777",
            "--domain-lb", "0",
            "--domain-ub", "3",
            "--treewidth-threshold", "3",
            "--output", output_path,
        ])

        @test result.config.count == 2
        @test result.config.nvars == [5]
        @test result.config.exact == false
        @test [row.strategy for row in result.rows] == [
            "primes",
            "prime_powers",
            "full",
        ]
        @test all(row -> row.constraints == 2, result.rows)
        @test all(row -> row.bounds_considered == 4, result.rows)
        @test all(row -> 0.0 <= row.fraction_bounds_tightened <= 1.0, result.rows)
        @test all(row -> 0.0 <= row.avg_relative_bound_reduction <= 1.0, result.rows)
        @test all(row -> row.total_residue_time_sec >= 0.0, result.rows)
        @test all(row -> row.avg_wall_time_sec_per_constraint >= 0.0, result.rows)
        @test isfile(output_path)

        header, rows = residue_family_csv_table(output_path)
        @test header == [
            "nvars",
            "modulus_limit",
            "strategy",
            "moduli",
            "num_moduli",
            "constraints",
            "bounds_considered",
            "bounds_tightened",
            "fraction_bounds_tightened",
            "avg_relative_bound_reduction",
            "total_residue_time_sec",
            "avg_wall_time_sec_per_constraint",
        ]
        @test length(rows) == 3
        @test all(row -> length(row) == length(header), rows)
    end
end

@testset "residue moduli family experiment writes exact summaries for all requested sizes" begin
    mktempdir() do dir
        output_path = joinpath(dir, "residue_family_summary_exact.csv")

        result = ResidueFamilyExperimentScript.main([
            "--count", "1",
            "--nvars", "4,5",
            "--seed-base", "888",
            "--domain-lb", "0",
            "--domain-ub", "1",
            "--treewidth-threshold", "3",
            "--exact=true",
            "--output", output_path,
        ])

        @test length(result.rows) == 6
        n4_rows = filter(row -> row.nvars == 4, result.rows)
        n5_rows = filter(row -> row.nvars == 5, result.rows)
        @test length(n4_rows) == 3
        @test length(n5_rows) == 3
        @test all(row -> row.exact_assignments_per_constraint == 16, n4_rows)
        @test all(row -> row.exact_assignments_per_constraint == 32, n5_rows)
        @test all(row -> 0.0 <= row.optimal_avg_relative_bound_reduction <= 1.0, n4_rows)
        @test all(row -> 0.0 <= row.fraction_bounds_tightened_to_optimal <= 1.0, n4_rows)
        @test all(row -> 0.0 <= row.optimal_avg_relative_bound_reduction <= 1.0, n5_rows)
        @test all(row -> 0.0 <= row.fraction_bounds_tightened_to_optimal <= 1.0, n5_rows)
        @test isfile(output_path)

        header, rows = residue_family_csv_table(output_path)
        @test header == [
            "nvars",
            "modulus_limit",
            "strategy",
            "moduli",
            "num_moduli",
            "constraints",
            "bounds_considered",
            "bounds_tightened",
            "fraction_bounds_tightened",
            "avg_relative_bound_reduction",
            "total_residue_time_sec",
            "avg_wall_time_sec_per_constraint",
            "exact_assignments_per_constraint",
            "optimal_avg_relative_bound_reduction",
            "bounds_tightened_to_optimal",
            "fraction_bounds_tightened_to_optimal",
        ]
        @test length(rows) == 6
        @test all(row -> length(row) == length(header), rows)
    end
end
