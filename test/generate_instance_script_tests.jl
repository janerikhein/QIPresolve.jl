using Test
import MathOptInterface as BenchmarkMOI
import QIPresolve as BenchmarkQIP

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
const INCREASING_QIP_TEST_INSTANCE_SCRIPT = joinpath(
    INSTANCE_SCRIPT_REPO_DIR,
    "scripts",
    "generate_increasing_random_qip_test_instances.jl",
)
const MAIN_BENCHMARK_INSTANCE_SCRIPT = joinpath(
    INSTANCE_SCRIPT_REPO_DIR,
    "scripts",
    "generate_main_benchmark_instances.jl",
)

include(MAIN_BENCHMARK_INSTANCE_SCRIPT)
const MainBenchmarkGenerator = Main.MainBenchmarkInstanceGenerator

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
    "box_margin",
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

const EXPECTED_INCREASING_QIP_TEST_INSTANCE_HEADER = [
    "instance_name",
    "file_name",
    "num",
    "created_at",
    "nvars",
    "ncons",
    "seed",
    "eq_constraints",
    "ineq_constraints",
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

@testset "main benchmark generator defines the complete specification" begin
    random = MainBenchmarkGenerator.random_specs()
    embedding = MainBenchmarkGenerator.embedding_specs()

    @test length(random) == 120
    @test length(embedding) == 150
    @test allunique(MainBenchmarkGenerator.random_instance_name.(random))
    @test allunique(MainBenchmarkGenerator.embedding_instance_name.(embedding))

    expected_random_probabilities = Dict(
        "bilinear" => (bilin = 0.2, diag = 0.0, lin = 0.0),
        "separable" => (bilin = 0.0, diag = 0.2, lin = 0.2),
        "pure" => (bilin = 0.2, diag = 0.2, lin = 0.0),
        "generic" => (bilin = 0.2, diag = 0.2, lin = 0.2),
    )
    for type in MainBenchmarkGenerator.RANDOM_TYPES
        family = filter(spec -> spec.type == type, random)
        @test length(family) == 30
        @test getproperty.(family, :seed) == collect(10_000:10_029)
        @test all(spec -> spec.probabilities == expected_random_probabilities[type], family)
    end

    for exactness in MainBenchmarkGenerator.EMBEDDING_EXACTNESSES
        for graph_type in MainBenchmarkGenerator.EMBEDDING_GRAPH_TYPES
            seed_base = MainBenchmarkGenerator.embedding_seed_base(exactness, graph_type)
            for anchoring in MainBenchmarkGenerator.EMBEDDING_ANCHORINGS
                family = filter(
                    spec -> spec.exactness == exactness &&
                        spec.graph_type == graph_type &&
                        spec.anchoring == anchoring,
                    embedding,
                )
                @test length(family) == 15
                @test getproperty.(family, :seed) == collect(seed_base:(seed_base + 14))
                @test all(spec -> spec.n == (exactness == "exact" ? 26 : 28), family)
                @test all(spec -> spec.alpha == (exactness == "exact" ? 0.0 : 0.01), family)
                @test all(
                    spec -> spec.num_anchors == (anchoring == "anchored" ? 3 : 0),
                    family,
                )
            end
        end
    end

    infeasible = filter(spec -> spec.feasibility == "infeasible", embedding)
    @test length(infeasible) == 30
    for strategy in MainBenchmarkGenerator.INFEASIBLE_EMBEDDING_STRATEGIES
        family = filter(spec -> spec.infeas_strategy == strategy, infeasible)
        @test length(family) == 15
        @test getproperty.(family, :seed) == collect(40_000:40_014)
        @test all(spec -> spec.exactness == "exact", family)
        @test all(spec -> spec.graph_type == "globally_rigid", family)
        @test all(spec -> spec.anchoring == "anchored", family)
        @test all(spec -> spec.n == 26, family)
        @test all(spec -> spec.num_anchors == 3, family)
        @test all(spec -> spec.alpha == 0.0, family)
    end
end

@testset "main benchmark generator writes family CSVs and protects its target" begin
    mktempdir() do dir
        target = joinpath(dir, "main_benchmark")
        random = [first(MainBenchmarkGenerator.random_specs())]
        all_embedding = MainBenchmarkGenerator.embedding_specs()
        feasible = first(filter(spec -> spec.feasibility == "feasible", all_embedding))
        bounding_box = first(filter(
            spec -> spec.infeas_strategy == :bounding_box,
            all_embedding,
        ))
        contraction = first(filter(
            spec -> spec.infeas_strategy == :vertex_contraction,
            all_embedding,
        ))
        embedding = [feasible, bounding_box, contraction]

        redirect_stdout(devnull) do
            MainBenchmarkGenerator.generate_dataset(
                target;
                random = random,
                embedding = embedding,
            )
        end

        random_header, random_rows = read_csv_table(joinpath(target, "random_instances.csv"))
        embedding_header, embedding_rows = read_csv_table(joinpath(target, "embedding_instances.csv"))
        @test length(random_rows) == 1
        @test length(embedding_rows) == 3
        @test csv_field(random_header, random_rows[1], "nvars") == "50"
        @test csv_field(random_header, random_rows[1], "ncons_requested") == "100"
        @test csv_field(random_header, random_rows[1], "p_con_eq") == "0.5"
        @test csv_field(random_header, random_rows[1], "var_threshold_lb") == "-10"
        @test csv_field(random_header, random_rows[1], "var_threshold_ub") == "10"
        @test csv_field(random_header, random_rows[1], "p_var_is_candidate") == "0.2"
        @test csv_field(random_header, random_rows[1], "force_bilin_even") == "true"
        @test csv_field(random_header, random_rows[1], "constraint_slack_range") == "-5:5"

        @test csv_field(embedding_header, embedding_rows[1], "R") == "100"
        @test csv_field(embedding_header, embedding_rows[1], "edge_density") == "0.05"
        @test csv_field(embedding_header, embedding_rows[1], "num_anchors") == "0"
        @test csv_field(embedding_header, embedding_rows[1], "alpha") == "0.0"
        @test csv_field(embedding_header, embedding_rows[1], "feasibility") == "feasible"
        @test csv_field(embedding_header, embedding_rows[2], "graph_type") == "globally_rigid"
        @test csv_field(embedding_header, embedding_rows[2], "num_anchors") == "3"
        @test csv_field(embedding_header, embedding_rows[2], "alpha") == "0.0"
        @test csv_field(embedding_header, embedding_rows[2], "feasibility") == "infeasible"
        @test csv_field(embedding_header, embedding_rows[2], "infeas_strategy") == "bounding_box"
        @test csv_field(embedding_header, embedding_rows[2], "infeas_base") == "globally_rigid"
        @test csv_field(embedding_header, embedding_rows[2], "box_margin") == "1"
        @test csv_field(embedding_header, embedding_rows[3], "infeas_strategy") == "vertex_contraction"
        @test csv_field(embedding_header, embedding_rows[3], "contraction_vertices") == "auto"
        @test csv_field(embedding_header, embedding_rows[2], "seed") == "40003"
        @test csv_field(embedding_header, embedding_rows[3], "seed") == "40000"

        lp_files = sort(filter(endswith(".lp"), readdir(target)))
        @test length(lp_files) == 4
        for file_name in lp_files
            model = BenchmarkMOI.FileFormats.Model(format = BenchmarkMOI.FileFormats.FORMAT_LP)
            BenchmarkMOI.read_from_file(model, joinpath(target, file_name))
            core_model = BenchmarkQIP.build_model(BenchmarkQIP.from_moi(model))
            @test all(var.lb <= var.ub for var in values(core_model.vars))
            if startswith(file_name, "random_")
                @test all(BenchmarkQIP.PresolvingCore.is_integer, core_model.cons)
            end
        end

        @test_throws ErrorException redirect_stdout(devnull) do
            MainBenchmarkGenerator.generate_dataset(
                target;
                random = random,
                embedding = embedding,
            )
        end


        write(joinpath(target, "obsolete"), "old benchmark")
        redirect_stdout(devnull) do
            MainBenchmarkGenerator.generate_dataset(
                target;
                random = random,
                embedding = embedding,
                force = true,
            )
        end
        @test !isfile(joinpath(target, "obsolete"))
        @test !ispath("$target.backup")
    end
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

@testset "graph embedding script retries rejected bounding-box seeds" begin
    mktempdir() do dir
        target = joinpath(dir, "graph_embedding_instances")
        csv_path = joinpath(target, "instances.csv")
        args = [
            "--type", "infeas",
            "--count", "1",
            "--target", target,
            "--csv", csv_path,
            "--n", "26",
            "--R", "100",
            "--num-anchors", "3",
            "--seed-base", "40000",
            "--infeas-strategy", "bounding_box",
            "--box-margin", "1",
        ]

        run_instance_script(GRAPH_INSTANCE_SCRIPT, args)
        run_instance_script(GRAPH_INSTANCE_SCRIPT, args)

        header, rows = read_csv_table(csv_path)
        @test String.(header) == EXPECTED_GRAPH_INSTANCE_HEADER
        @test csv_field(header, rows[1], "seed") == "40003"
        @test csv_field(header, rows[2], "seed") == "40004"
        @test all(row -> csv_field(header, row, "box_margin") == "1", rows)

        legacy_args = copy(args)
        legacy_args[end - 1] = "--box-scale"
        @test_throws Base.ProcessFailedException run_instance_script(
            GRAPH_INSTANCE_SCRIPT, legacy_args
        )
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

@testset "increasing random QIP test instance script writes 10 mixed-constraint instances" begin
    mktempdir() do dir
        target = joinpath(dir, "test_instances")
        csv_path = joinpath(target, "instances.csv")
        seed_base = 3000
        max_attempts = 100
        args = [
            "--target", target,
            "--csv", csv_path,
            "--seed-base", string(seed_base),
            "--max-attempts", string(max_attempts),
        ]

        run_instance_script(INCREASING_QIP_TEST_INSTANCE_SCRIPT, args)

        header, rows = read_csv_table(csv_path)
        @test String.(header) == EXPECTED_INCREASING_QIP_TEST_INSTANCE_HEADER
        @test length(rows) == 10
        @test all(row -> length(row) == length(header), rows)

        for (idx, row) in enumerate(rows)
            nvars = 10 * idx
            ncons = 2 * nvars
            expected_name = "test_qip_n$nvars"
            expected_file = "$expected_name.lp"
            seed = parse(Int, csv_field(header, row, "seed"))
            eq_constraints = parse(Int, csv_field(header, row, "eq_constraints"))
            ineq_constraints = parse(Int, csv_field(header, row, "ineq_constraints"))

            @test csv_field(header, row, "instance_name") == expected_name
            @test csv_field(header, row, "file_name") == expected_file
            @test csv_field(header, row, "num") == string(idx)
            @test csv_field(header, row, "nvars") == string(nvars)
            @test csv_field(header, row, "ncons") == string(ncons)
            @test seed_base + (idx - 1) * max_attempts <= seed
            @test seed < seed_base + idx * max_attempts
            @test eq_constraints > 0
            @test ineq_constraints > 0
            @test eq_constraints + ineq_constraints == ncons
            @test csv_field(header, row, "p_con_eq") == "0.5"
            @test isfile(joinpath(target, expected_file))
        end
    end
end
