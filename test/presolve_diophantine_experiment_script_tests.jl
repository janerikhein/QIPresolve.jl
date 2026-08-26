using Random
using Test
import QIPresolve.PresolvingCore as PC

include(joinpath(@__DIR__, "..", "scripts", "presolve_diophantine_experiment.jl"))
const DiophantineExperimentScript = Main.DiophantineParityPresolveExperiment

const DiophantineExperimentQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const DiophantineExperimentLinTerm = Tuple{Float64, PC.VarId}

function diophantine_experiment_empty_objective()
    return PC.QuadExpr(DiophantineExperimentQuadTerm[], DiophantineExperimentLinTerm[])
end

function diophantine_experiment_term_flags(model::PC.QPModel)
    has_linear = false
    has_diagonal = false
    has_bilinear = false

    for con in model.cons
        var_ids = collect(PC.vars(con.qe))
        for var_id in var_ids
            isapprox(PC.get_lin_coeff(con.qe, var_id), 0.0) || (has_linear = true)
            isapprox(PC.get_quad_coeff(con.qe, var_id, var_id), 0.0) || (has_diagonal = true)
        end

        for (index, first_id) in enumerate(var_ids)
            for second_index in (index + 1):length(var_ids)
                second_id = var_ids[second_index]
                isapprox(PC.get_quad_coeff(con.qe, first_id, second_id), 0.0) ||
                    (has_bilinear = true)
            end
        end
    end

    return (
        linear = has_linear,
        diagonal = has_diagonal,
        bilinear = has_bilinear,
    )
end

@testset "diophantine experiment parses type probabilities" begin
    @test DiophantineExperimentScript.canonical_type("seperable") == "separable"
    @test DiophantineExperimentScript.canonical_type("purely_quadratic") == "pure"
    @test DiophantineExperimentScript.parse_types("bilinear,seperable,general") ==
          ["bilinear", "separable", "general"]
    @test DiophantineExperimentScript.parse_parity_strategy("mod2-basic") == :mod2_basic
    @test DiophantineExperimentScript.parse_parity_strategy("mod4_basic") == :mod4_basic
    @test DiophantineExperimentScript.parse_parity_strategy("full") == :full
    @test_throws ErrorException DiophantineExperimentScript.parse_parity_strategy("basic")

    p = 0.25
    @test DiophantineExperimentScript.type_probabilities("bilinear", p) ==
          (p_var_bilin = p, p_var_diag = 0.0, p_var_lin = 0.0)
    @test DiophantineExperimentScript.type_probabilities("separable", p) ==
          (p_var_bilin = 0.0, p_var_diag = p, p_var_lin = p)
    @test DiophantineExperimentScript.type_probabilities("pure", p) ==
          (p_var_bilin = p, p_var_diag = p, p_var_lin = 0.0)
    @test DiophantineExperimentScript.type_probabilities("general", p) ==
          (p_var_bilin = p, p_var_diag = p, p_var_lin = p)
end

@testset "diophantine experiment generator respects class support" begin
    expected_flags = Dict(
        "bilinear" => (linear = false, diagonal = false, bilinear = true),
        "separable" => (linear = true, diagonal = true, bilinear = false),
        "pure" => (linear = false, diagonal = true, bilinear = true),
        "general" => (linear = true, diagonal = true, bilinear = true),
    )

    for (index, type) in enumerate(DiophantineExperimentScript.DEFAULT_TYPES)
        model, upper_bounds, x_star = DiophantineExperimentScript.build_random_diophantine_model(
            type,
            4,
            2,
            1.0,
            1000 + index,
            MersenneTwister(2000 + index),
        )

        @test diophantine_experiment_term_flags(model) == expected_flags[type]
        @test collect(PC.vars(model.obj_expr)) == PC.VarId[]
        @test all(model.vars[var_id].lb == 0.0 for var_id in 1:4)
        @test [Int(model.vars[var_id].ub) for var_id in 1:4] == upper_bounds
        @test all(con -> con.lhs == con.rhs, model.cons)
        @test all(con -> PC.eval_full(con.qe, x_star) == con.rhs, model.cons)
    end
end

@testset "diophantine experiment recursively classifies postsolve bits" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        3 => PC.IntVar(0.0, 1.0),
    )
    model = PC.QPModel(vars, PC.Constraint[], diophantine_experiment_empty_objective(), :min)
    postsolver = PC.ParityPostsolver([1])
    PC.append_binary_bit!(postsolver, 1, 3)

    pattern_counts = DiophantineExperimentScript.postsolve_bit_counts(postsolver, model, [3])
    @test pattern_counts.fixed_bits == 0
    @test pattern_counts.pattern_bits == 1
    @test pattern_counts.total_bits == 2

    PC.register_fixed_var!(postsolver, 3, 1.0)
    fixed_ref_counts = DiophantineExperimentScript.postsolve_bit_counts(postsolver, model, [3])
    @test fixed_ref_counts.fixed_bits == 1
    @test fixed_ref_counts.pattern_bits == 0
    @test fixed_ref_counts.total_bits == 2

    PC.register_fixed_var!(postsolver, 1, 0.0)
    fixed_variable_counts = DiophantineExperimentScript.postsolve_bit_counts(postsolver, model, [3])
    @test fixed_variable_counts.fixed_bits == 2
    @test fixed_variable_counts.pattern_bits == 0
    @test fixed_variable_counts.total_bits == 2
end

@testset "diophantine experiment dense parity state comparison is stack-safe" begin
    vars = Dict{PC.VarId, PC.IntVar}(var_id => PC.IntVar(0.0, 100.0) for var_id in 1:100)
    quad_terms = DiophantineExperimentQuadTerm[]
    for first_id in 1:100
        for second_id in (first_id + 1):100
            push!(quad_terms, (1.0, first_id, second_id))
            length(quad_terms) >= 1500 && break
        end
        length(quad_terms) >= 1500 && break
    end

    cons = [
        PC.Constraint(
            con_id,
            PC.QuadExpr(copy(quad_terms), DiophantineExperimentLinTerm[]),
            0.0,
            0.0,
        )
        for con_id in 1:25
    ]
    model = PC.QPModel(vars, cons, diophantine_experiment_empty_objective(), :min)

    @test PC._model_state_signature(model) == PC._model_state_signature(model)

    same_cons = [
        PC.Constraint(
            con_id,
            PC.QuadExpr(copy(quad_terms), DiophantineExperimentLinTerm[]),
            Float64(con_id - 1),
            Float64(30 - con_id),
        )
        for con_id in 1:2
    ]
    same_model = PC.QPModel(
        deepcopy(vars),
        same_cons,
        diophantine_experiment_empty_objective(),
        :min,
    )
    PC.normalize!(same_model)
    @test !same_model.infeasible
    @test length(same_model.cons) == 1
end

@testset "diophantine experiment writes summary and tuple outputs" begin
    mktempdir() do dir
        stdout_path = joinpath(dir, "stdout.txt")
        result = open(stdout_path, "w") do io
            redirect_stdout(io) do
                DiophantineExperimentScript.main([
                    "--output-dir", dir,
                    "--nvars", "4",
                    "--instances-per-cell", "1",
                    "--constraints", "1",
                    "--types", "bilinear,seperable",
                    "--seed-base", "41000",
                    "--parity-strategy", "mod2-basic",
                ])
            end
        end

        @test result.config.parity_strategy == :mod2_basic
        @test result.instance_count == 2
        @test length(result.rows) == 2
        @test result.summary_path == joinpath(dir, "parity_presolve_summary.csv")

        printed_lines = readlines(stdout_path)
        @test length(printed_lines) == 2
        @test startswith(printed_lines[1], "bilinear,1,")
        @test startswith(printed_lines[2], "separable,1,")

        summary_lines = readlines(result.summary_path)
        @test first(summary_lines) == DiophantineExperimentScript.SUMMARY_HEADER
        @test length(summary_lines) == 3
        @test summary_lines[2:end] == printed_lines

        for type in ["bilinear", "separable"]
            density_lines = readlines(joinpath(dir, "density_vs_dom_red_$type.txt"))
            domain_lines = readlines(joinpath(dir, "avg_domain_size_vs_dom_red_$type.txt"))
            @test length(density_lines) == 1
            @test length(domain_lines) == 1
            @test startswith(first(density_lines), "(")
            @test startswith(first(domain_lines), "(")
        end
    end
end
