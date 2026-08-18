using Test
using JuMP: backend
import JSON
import QIPresolve.PresolvingCore as PC
using QIPresolve.InstanceGeneration: generate_random_qip_model
using QIPresolve.ModelIO: save_moi

include(joinpath(@__DIR__, "..", "scripts", "presolve_lp_scip_stats.jl"))
const LPStatsScript = Main.PresolveLpScipStatsScript

@testset "presolve LP SCIP stats script writes JSON" begin
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

        result_path = LPStatsScript.main([
            lp_path,
            "--output",
            output_path,
            "--silent",
        ])

        @test result_path == output_path
        @test isfile(output_path)

        data = JSON.parsefile(output_path)
        required_keys = [
            "instance_name",
            "metadata",
            "parameters",
            "scip_original",
            "scip_presolved",
            "model_loaded",
            "model_before",
            "model_after",
            "parity_stats",
            "residue_stats",
            "validation",
        ]
        for key in required_keys
            @test haskey(data, key)
        end

        @test data["instance_name"] == "tiny"
        @test data["metadata"]["legacy_lp_repaired"] == false
        @test data["metadata"]["legacy_lp_repaired_constraints"] == 0
        @test Set(keys(data["scip_original"])) == Set(keys(data["scip_presolved"]))
        @test haskey(data["scip_original"], "termination_status")
        @test haskey(data["scip_original"], "node_count")
        @test haskey(data["scip_presolved"], "termination_status")
        @test haskey(data["scip_presolved"], "node_count")

        @test haskey(data["model_before"], "num_variables")
        @test haskey(data["model_after"], "num_variables")
        @test haskey(data["parity_stats"], "num_parity_constraints_generated")
        @test haskey(data["residue_stats"], "num_constraints_processed")

        @test data["validation"]["checked"] == true
        @test data["validation"]["feasible"] == true
        @test data["validation"]["constraint_violation_count"] == 0
    end
end

@testset "presolve LP SCIP stats script repairs legacy ranged quadratic LP" begin
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

        result_path = LPStatsScript.main([
            lp_path,
            "--output",
            output_path,
            "--silent",
        ])

        @test result_path == output_path
        data = JSON.parsefile(output_path)
        @test data["metadata"]["legacy_lp_repaired"] == true
        @test data["metadata"]["legacy_lp_repaired_constraints"] == 1
        @test data["metadata"]["legacy_lp_repair_error"] isa String
        @test data["model_loaded"]["num_constraints"] == 2
    end
end

@testset "random QIP LP generation writes reader-compatible LP" begin
    mktempdir() do dir
        lp_path = joinpath(dir, "generated_random_qip.lp")
        jump_model, _ = generate_random_qip_model(
            4,
            2;
            p_con_eq = 0.0,
            var_threshold_lb = -3,
            var_threshold_ub = 3,
            p_var_is_candidate = 1.0,
            p_var_bilin = 0.5,
            p_var_diag = 0.5,
            p_var_lin = 0.0,
            coeff_lb = -5,
            coeff_ub = 5,
            force_diag_even = false,
            force_lin_even = false,
            force_feasibility = true,
            constraint_slack_range = [-1, 1],
            seed = 17,
        )

        save_moi(backend(jump_model), lp_path)
        contents = read(lp_path, String)

        @test !occursin(r":\s*-?[\d.]+(?:[eE][+-]?\d+)?\s*<=\s*\[", contents)

        load_info = LPStatsScript.load_lp_core_model(lp_path)
        @test load_info.model isa PC.QPModel
        @test load_info.legacy_lp_repaired == false
        @test load_info.legacy_lp_repaired_constraints == 0
    end
end
