using Test
import JSON
import QIPresolve.PresolvingCore as PC

include(joinpath(@__DIR__, "..", "scripts", "presolve_lp_stats.jl"))
const LPOnlyStatsScript = Main.PresolveLpStatsScript

@testset "presolve LP stats script writes presolve-only JSON" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "tiny.lp")
        output_path = joinpath(dir, "tiny_stats.json")

        write(
            lp_path,
            """
            Maximize
             obj: x1 + x2
            Subject To
             c1: x1 + x2 <= 1
            Binary
             x1
             x2
            End
            """,
        )

        result_path = LPOnlyStatsScript.main([
            lp_path,
            "--output",
            output_path,
        ])

        @test result_path == output_path
        @test isfile(output_path)

        data = JSON.parsefile(output_path)
        required_keys = [
            "instance_name",
            "metadata",
            "parameters",
            "model_loaded",
            "model_before",
            "model_after",
            "parity_stats",
            "residue_stats",
        ]
        for key in required_keys
            @test haskey(data, key)
        end

        omitted_keys = [
            "scip_original",
            "scip_presolved",
            "validation",
        ]
        for key in omitted_keys
            @test !haskey(data, key)
        end

        @test data["instance_name"] == "tiny"
        @test data["metadata"]["legacy_lp_repaired"] == false
        @test data["metadata"]["legacy_lp_repaired_constraints"] == 0
        @test data["metadata"]["script"] == "presolve_lp_stats.jl"
        @test !haskey(data["metadata"], "scip_version")

        @test haskey(data["parameters"], "residue_strategy")
        @test haskey(data["parameters"], "residue_threshold")
        @test haskey(data["parameters"], "treewidth_threshold")
        @test !haskey(data["parameters"], "scip_config")
        @test !haskey(data["parameters"], "silent")
        @test !haskey(data["parameters"], "validation_tolerance")

        @test haskey(data["model_before"], "num_variables")
        @test haskey(data["model_after"], "num_variables")
        @test haskey(data["parity_stats"], "num_parity_constraints_generated")
        @test haskey(data["residue_stats"], "num_constraints_processed")
    end
end

@testset "presolve LP stats script repairs legacy ranged quadratic LP" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "legacy_ranged.lp")
        output_path = joinpath(dir, "legacy_ranged_stats.json")

        write(
            lp_path,
            """
            Minimize
             obj: x1 + x2
            Subject To
             c1: 0 <= [ 1 x1 ^ 2 + 2 x1 * x2 ] <= 1
            Binary
             x1
             x2
            End
            """,
        )

        result_path = LPOnlyStatsScript.main([
            lp_path,
            "--output",
            output_path,
        ])

        @test result_path == output_path
        data = JSON.parsefile(output_path)
        @test data["metadata"]["legacy_lp_repaired"] == true
        @test data["metadata"]["legacy_lp_repaired_constraints"] == 1
        @test data["metadata"]["legacy_lp_repair_error"] isa String
        @test data["model_loaded"]["num_constraints"] == 2
    end
end

@testset "presolve LP stats script rejects SCIP-only options" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "tiny.lp")

        write(
            lp_path,
            """
            Maximize
             obj: x1
            Binary
             x1
            End
            """,
        )

        @test_throws ErrorException LPOnlyStatsScript.parse_args([
            lp_path,
            "--scip-config",
            joinpath(dir, "settings.set"),
        ])
    end
end
