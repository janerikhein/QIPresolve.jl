using Test
import QIPresolve.PresolvingCore as PC

include(joinpath(@__DIR__, "..", "scripts", "presolve_lattice_experiment.jl"))
const LatticeExperimentScript = Main.LatticeParityPresolveExperiment

const LatticeExperimentQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LatticeExperimentLinTerm = Tuple{Float64, PC.VarId}

function lattice_experiment_empty_objective()
    return PC.QuadExpr(LatticeExperimentQuadTerm[], LatticeExperimentLinTerm[])
end

@testset "lattice experiment parses CLI and type aliases" begin
    @test LatticeExperimentScript.canonical_type("sparse") == "2-connected-sparse"
    @test LatticeExperimentScript.canonical_type("con_dense") == "2-connected-dense"
    @test LatticeExperimentScript.canonical_type("lam") == "laman"
    @test LatticeExperimentScript.canonical_type("gr") == "globally-rigid"
    @test LatticeExperimentScript.parse_types("sparse,dense,lam,gr,sparse") == [
        "2-connected-sparse",
        "2-connected-dense",
        "laman",
        "globally-rigid",
    ]
    @test LatticeExperimentScript.parse_parity_strategy("mod2-basic") == :mod2_basic
    @test LatticeExperimentScript.parse_parity_strategy("mod4_basic") == :mod4_basic
    @test LatticeExperimentScript.parse_parity_strategy("full") == :full
    @test_throws ErrorException LatticeExperimentScript.parse_parity_strategy("basic")

    config = LatticeExperimentScript.build_config([
        "--num-anchors", "3",
        "--instances", "2",
        "--types", "sparse,gr",
        "--seed-base", "7",
        "--seed-step", "2",
        "--output-dir", "tmp-lattice-results",
        "--parity-strategy", "mod4-basic",
        "--free-vertices", "5",
    ])
    @test config.num_anchors == 3
    @test config.instances_per_type == 2
    @test config.types == ["2-connected-sparse", "globally-rigid"]
    @test config.seed_base == 7
    @test config.seed_step == 2
    @test config.parity_strategy == :mod4_basic
    @test config.free_vertices == 5
    @test LatticeExperimentScript.effective_n(config) == 8
end

@testset "lattice experiment computes effective n and graph edge counts" begin
    @test LatticeExperimentScript.effective_n(0) == 51
    @test LatticeExperimentScript.effective_n(2) == 52
    @test LatticeExperimentScript.effective_n(0, 3) == 4
    @test LatticeExperimentScript.effective_n(2, 3) == 5

    n = 51
    m_min = n
    m_max = n * (n - 1) ÷ 2
    @test LatticeExperimentScript.lattice_edge_count("2-connected-sparse", n) ==
          round(Int, m_min + 0.05 * (m_max - m_min))
    @test LatticeExperimentScript.lattice_edge_count("2-connected-dense", n) ==
          round(Int, m_min + 0.3 * (m_max - m_min))
    @test LatticeExperimentScript.lattice_edge_count("laman", n) == 2n - 3
    @test LatticeExperimentScript.lattice_edge_count("globally-rigid", n) == 2n - 2
end

@testset "lattice experiment recursively classifies postsolve bits" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 3.0),
        3 => PC.IntVar(0.0, 1.0),
    )
    model = PC.QPModel(vars, PC.Constraint[], lattice_experiment_empty_objective(), :min)
    postsolver = PC.ParityPostsolver([1])
    bit_widths = Dict{PC.VarId, Int}(1 => 2)
    PC.append_binary_bit!(postsolver, 1, 3)

    pattern_counts = LatticeExperimentScript.postsolve_bit_counts(postsolver, model, bit_widths)
    @test pattern_counts.fixed_bits == 0
    @test pattern_counts.pattern_bits == 1
    @test pattern_counts.total_bits == 2

    PC.register_fixed_var!(postsolver, 3, 1.0)
    fixed_ref_counts = LatticeExperimentScript.postsolve_bit_counts(postsolver, model, bit_widths)
    @test fixed_ref_counts.fixed_bits == 1
    @test fixed_ref_counts.pattern_bits == 0
    @test fixed_ref_counts.total_bits == 2

    PC.register_fixed_var!(postsolver, 1, 0.0)
    fixed_variable_counts = LatticeExperimentScript.postsolve_bit_counts(postsolver, model, bit_widths)
    @test fixed_variable_counts.fixed_bits == 2
    @test fixed_variable_counts.pattern_bits == 0
    @test fixed_variable_counts.total_bits == 2
end

@testset "lattice experiment writes summary without instance files" begin
    mktempdir() do dir
        stdout_path = joinpath(dir, "stdout.txt")
        result = open(stdout_path, "w") do io
            redirect_stdout(io) do
                LatticeExperimentScript.main([
                    "--output-dir", dir,
                    "--instances-per-type", "1",
                    "--types", "sparse,laman",
                    "--seed-base", "51000",
                    "--parity-strategy", "mod2-basic",
                    "--free-vertices", "3",
                ])
            end
        end

        @test result.config.parity_strategy == :mod2_basic
        @test result.instance_count == 2
        @test length(result.rows) == 2
        @test result.summary_path == joinpath(dir, "lattice_parity_presolve_summary.csv")

        printed_lines = readlines(stdout_path)
        @test length(printed_lines) == 2
        @test startswith(printed_lines[1], "2-connected-sparse,")
        @test startswith(printed_lines[2], "laman,")

        summary_lines = readlines(result.summary_path)
        @test first(summary_lines) == LatticeExperimentScript.SUMMARY_HEADER
        @test length(summary_lines) == 3
        @test summary_lines[2:end] == printed_lines
        @test isempty(filter(path -> endswith(path, ".lp"), readdir(dir; join = true)))
    end
end
