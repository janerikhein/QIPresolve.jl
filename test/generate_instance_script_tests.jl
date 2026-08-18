using Test

const INSTANCE_SCRIPT_REPO_DIR = normpath(joinpath(@__DIR__, ".."))
const GRAPH_INSTANCE_SCRIPT = joinpath(
    INSTANCE_SCRIPT_REPO_DIR,
    "scripts",
    "generate_graph_embedding_instances.jl",
)
const QIP_INSTANCE_SCRIPT = joinpath(
    INSTANCE_SCRIPT_REPO_DIR,
    "scripts",
    "generate_random_qip_instances.jl",
)

const EXPECTED_GRAPH_INSTANCE_HEADER = [
    "instance_name",
    "type",
    "num",
    "created_at",
    "n",
    "R",
    "seed",
    "num_anchors",
    "alpha",
    "edge_density",
    "pH2",
    "max_coord_tries",
    "max_tries_H2",
    "infeas_strategy",
    "infeas_base",
    "box_scale",
]

const EXPECTED_QIP_INSTANCE_HEADER = [
    "instance_name",
    "num",
    "created_at",
    "nvars",
    "ncons",
    "seed",
    "p_con_eq",
    "var_threshold_lb",
    "var_threshold_ub",
    "p_var_is_candidate",
    "p_var_bilin",
    "p_var_diag",
    "p_var_lin",
    "coeff_lb",
    "coeff_ub",
    "force_diag_even",
    "force_lin_even",
    "force_feasibility",
    "constraint_slack_range",
]

function run_instance_script(script_path::AbstractString, args::Vector{String})
    cmd = `$(Base.julia_cmd()) --project=$INSTANCE_SCRIPT_REPO_DIR $script_path $args`
    run(pipeline(cmd; stdout = devnull))
    return nothing
end

split_csv_line(line::AbstractString) = split(line, ','; keepempty = true)

function read_csv_table(path::AbstractString)
    lines = readlines(path)
    return split_csv_line(first(lines)), split_csv_line.(lines[2:end])
end

function csv_field(header::Vector{SubString{String}}, row::Vector{SubString{String}}, name::String)
    idx = findfirst(==(name), header)
    idx === nothing && error("missing CSV column $name")
    return row[idx]
end

@testset "graph embedding instance script writes slim CSV and continues numbering" begin
    mktempdir() do dir
        target = joinpath(dir, "graph_embedding_instances")
        csv_path = joinpath(target, "instances.csv")
        args = [
            "--type", "con",
            "--count", "1",
            "--target", target,
            "--csv", csv_path,
            "--instance-prefix", "sample_",
            "--n", "5",
            "--R", "8",
            "--edge-density", "0.25",
            "--seed-base", "1000",
        ]

        run_instance_script(GRAPH_INSTANCE_SCRIPT, args)
        run_instance_script(GRAPH_INSTANCE_SCRIPT, args)

        header, rows = read_csv_table(csv_path)
        @test String.(header) == EXPECTED_GRAPH_INSTANCE_HEADER
        @test length(rows) == 2
        @test all(row -> length(row) == length(header), rows)
        @test csv_field(header, rows[1], "instance_name") == "sample_con_1"
        @test csv_field(header, rows[2], "instance_name") == "sample_con_2"
        @test csv_field(header, rows[1], "num") == "1"
        @test csv_field(header, rows[2], "num") == "2"
        @test csv_field(header, rows[1], "seed") == "1000"
        @test csv_field(header, rows[2], "seed") == "1001"
        @test csv_field(header, rows[1], "edge_density") == "0.25"
        @test isfile(joinpath(target, "sample_con_1.lp"))
        @test isfile(joinpath(target, "sample_con_2.lp"))
    end
end

@testset "random QIP instance script writes slim CSV and continues numbering" begin
    mktempdir() do dir
        target = joinpath(dir, "random_qip_instances")
        csv_path = joinpath(target, "instances.csv")
        args = [
            "--count", "1",
            "--target", target,
            "--csv", csv_path,
            "--instance-prefix", "sample_",
            "--nvars", "4",
            "--ncons", "2",
            "--seed-base", "2000",
            "--p-var-is-candidate", "1.0",
            "--coeff-lb", "-5",
            "--coeff-ub", "5",
            "--constraint-slack-range", "-1:1",
        ]

        run_instance_script(QIP_INSTANCE_SCRIPT, args)
        run_instance_script(QIP_INSTANCE_SCRIPT, args)

        header, rows = read_csv_table(csv_path)
        @test String.(header) == EXPECTED_QIP_INSTANCE_HEADER
        @test length(rows) == 2
        @test all(row -> length(row) == length(header), rows)
        @test csv_field(header, rows[1], "instance_name") == "sample_qip_1"
        @test csv_field(header, rows[2], "instance_name") == "sample_qip_2"
        @test csv_field(header, rows[1], "num") == "1"
        @test csv_field(header, rows[2], "num") == "2"
        @test csv_field(header, rows[1], "seed") == "2000"
        @test csv_field(header, rows[2], "seed") == "2001"
        @test csv_field(header, rows[1], "constraint_slack_range") == "-1:1"
        @test isfile(joinpath(target, "sample_qip_1.lp"))
        @test isfile(joinpath(target, "sample_qip_2.lp"))
    end
end
