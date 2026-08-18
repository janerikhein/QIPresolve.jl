using Test
import CSV
import JSON

include(joinpath(@__DIR__, "..", "scripts", "merge_graph_embedding_presolve_stats.jl"))
const MergePresolveStatsScript = Main.MergeGraphEmbeddingPresolveStats

const EXPECTED_MERGED_PRESOLVE_HEADER = [
    "instance_name",
    "type",
    "n",
    "R",
    "alpha",
    "edge_density",
    "infeas_strategy",
    "infeas_base",
    "box_scale",
    "log_domain_sum_orig",
    "log_domain_sum_ps",
    "num_quadratic_constraints_orig",
    "num_quadratic_constraints_ps",
    "num_linear_constraints_orig",
    "num_linear_constraints_ps",
    "parity_presolve_time",
    "residue_presolve_time",
    "num_parity_substitutions",
    "num_pattern_substitutions",
    "num_variables_fixed",
    "num_constraints_tightened",
    "num_constraints_eq_tightened",
    "avg_relative_bound_tightening",
    "infeasibility_detected",
    "infeasibility_source",
]

function write_synthetic_presolve_stats(
        path::AbstractString,
        instance_name::AbstractString;
        parity_infeasible::Bool = false,
        parity_source::AbstractString = "none",
        residue_infeasible::Bool = false,
    )
    data = Dict{String, Any}(
        "instance_name" => instance_name,
        "model_loaded" => Dict{String, Any}(
            "log_domain_sum" => 10.5,
            "num_quadratic_constraints" => 7,
            "num_linear_constraints" => 3,
        ),
        "model_after" => Dict{String, Any}(
            "log_domain_sum" => 4.25,
            "num_quadratic_constraints" => 2,
            "num_linear_constraints" => 8,
        ),
        "parity_stats" => Dict{String, Any}(
            "parity_presolve_time" => 0.125,
            "num_parity_substitutions" => 11,
            "num_pattern_substitutions" => 5,
            "num_variables_fixed_after_parity_substitution" => 2,
            "infeasibility_detected" => parity_infeasible,
            "infeasibility_source" => parity_source,
        ),
        "residue_stats" => Dict{String, Any}(
            "residue_presolve_time" => 0.75,
            "num_constraints_tightened" => 13,
            "num_constraints_eq_tightened" => 1,
            "avg_relative_bound_tightening" => 0.375,
            "infeasibility_detected" => residue_infeasible,
        ),
    )

    open(path, "w") do io
        JSON.print(io, data)
        println(io)
    end
    return path
end

@testset "merge graph embedding presolve stats script writes requested CSV" begin
    mktempdir() do dir
        stats_dir = joinpath(dir, "presolve_stats")
        mkpath(stats_dir)
        instances_csv = joinpath(dir, "instances.csv")
        output_csv = joinpath(dir, "instances_with_presolve_stats.csv")

        write(
            instances_csv,
            """
            instance_name,type,num,created_at,n,R,seed,num_anchors,alpha,edge_density,pH2,max_coord_tries,max_tries_H2,infeas_strategy,infeas_base,box_scale
            residue_case,2_connected,1,2026-08-18T22:41:22,20,100,1,0,0.0,0.2,,10000,,,,
            parity_case,likely_infeasible,2,2026-08-18T22:41:22,30,100,2,0,0.0,,,,bounding_box,globally_rigid,0.75
            """,
        )

        write_synthetic_presolve_stats(
            joinpath(stats_dir, "residue_case_stats.json"),
            "residue_case";
            residue_infeasible = true,
        )
        write_synthetic_presolve_stats(
            joinpath(stats_dir, "parity_case_stats.json"),
            "parity_case";
            parity_infeasible = true,
            parity_source = "elimination",
            residue_infeasible = true,
        )

        result_path = MergePresolveStatsScript.main([
            "--instances-csv",
            instances_csv,
            "--stats-dir",
            stats_dir,
            "--output",
            output_csv,
        ])

        @test result_path == output_csv
        @test isfile(output_csv)

        output = CSV.File(output_csv)
        @test String.(propertynames(output)) == EXPECTED_MERGED_PRESOLVE_HEADER
        rows = collect(output)
        @test length(rows) == 2

        @test rows[1].instance_name == "residue_case"
        @test rows[1].type == "2_connected"
        @test rows[1].n == 20
        @test rows[1].R == 100
        @test rows[1].alpha == 0.0
        @test rows[1].edge_density == 0.2
        @test rows[1].log_domain_sum_orig == 10.5
        @test rows[1].log_domain_sum_ps == 4.25
        @test rows[1].num_quadratic_constraints_orig == 7
        @test rows[1].num_quadratic_constraints_ps == 2
        @test rows[1].num_linear_constraints_orig == 3
        @test rows[1].num_linear_constraints_ps == 8
        @test rows[1].parity_presolve_time == 0.125
        @test rows[1].residue_presolve_time == 0.75
        @test rows[1].num_parity_substitutions == 11
        @test rows[1].num_pattern_substitutions == 5
        @test rows[1].num_variables_fixed == 2
        @test rows[1].num_constraints_tightened == 13
        @test rows[1].num_constraints_eq_tightened == 1
        @test rows[1].avg_relative_bound_tightening == 0.375
        @test rows[1].infeasibility_detected == true
        @test rows[1].infeasibility_source == "residue"

        @test rows[2].infeasibility_detected == true
        @test rows[2].infeasibility_source == "elimination"
    end
end

@testset "merge graph embedding presolve stats script requires each stats file" begin
    mktempdir() do dir
        stats_dir = joinpath(dir, "presolve_stats")
        mkpath(stats_dir)
        instances_csv = joinpath(dir, "instances.csv")

        write(
            instances_csv,
            """
            instance_name,type,num,created_at,n,R,seed,num_anchors,alpha,edge_density,pH2,max_coord_tries,max_tries_H2,infeas_strategy,infeas_base,box_scale
            missing_case,2_connected,1,2026-08-18T22:41:22,20,100,1,0,0.0,0.2,,10000,,,,
            """,
        )

        @test_throws ErrorException MergePresolveStatsScript.main([
            "--instances-csv",
            instances_csv,
            "--stats-dir",
            stats_dir,
            "--output",
            joinpath(dir, "output.csv"),
        ])
    end
end
