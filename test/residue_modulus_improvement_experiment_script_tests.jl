using Test

const RESIDUE_MODULUS_EXPERIMENT_SCRIPT = joinpath(
    @__DIR__,
    "..",
    "scripts",
    "residue_modulus_improvement_experiment.jl",
)

include(RESIDUE_MODULUS_EXPERIMENT_SCRIPT)
const ResidueModulusExperimentScript = Main.ResidueModulusImprovementExperiment

function residue_experiment_csv_table(path::AbstractString)
    lines = readlines(path)
    return split(first(lines), ','), split.(lines[2:end], ','; keepempty = true)
end

@testset "residue modulus improvement experiment writes strategy summaries" begin
    mktempdir() do dir
        output_path = joinpath(dir, "residue_strategy_summary.csv")

        result = ResidueModulusExperimentScript.main([
            "--count", "1",
            "--nvars", "8",
            "--ncons", "4",
            "--seed-base", "321",
            "--max-modulus", "4",
            "--treewidth-threshold", "2",
            "--output", output_path,
        ])

        @test result.config.max_modulus == 4
        @test result.config.p_var_is_candidate == 0.1
        @test result.config.p_var_bilin == 0.3
        @test result.config.p_var_diag == 0.3
        @test result.config.p_var_lin == 0.3
        @test [row.strategy for row in result.rows] == [
            "all_2_M",
            "primes_2_M",
            "prime_powers_2_M",
        ]
        @test [row.num_moduli for row in result.rows] == [3, 2, 3]
        @test all(row -> row.strategy_time_sec >= 0.0, result.rows)
        @test length(result.rows) == 3
        @test length(result.csv_rows) == 8
        @test isfile(output_path)

        header, rows = residue_experiment_csv_table(output_path)
        @test header == [
            "strategy",
            "max_modulus",
            "modulus_order",
            "modulus",
            "num_moduli",
            "processed_constraints",
            "improved_constraints",
            "avg_relative_bound_tightening",
            "strategy_time_sec",
            "infeasibilities",
            "modulus_evaluated_constraints",
            "modulus_bound_tightened_constraints",
            "modulus_lower_bound_tightenings",
            "modulus_upper_bound_tightenings",
            "modulus_infeasibilities",
        ]
        @test length(rows) == 8
        @test all(row -> length(row) == length(header), rows)
    end
end
