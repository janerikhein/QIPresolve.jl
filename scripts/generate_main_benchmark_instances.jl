#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
end

module MainBenchmarkInstanceGenerator

using CSV
import JuMP
import MathOptInterface as MOI

using QIPresolve.InstanceGeneration:
    BoundingBoxCandidateRejected,
    generate_2_connected_instance,
    generate_laman_instance,
    generate_likely_infeasible_embedding_instance,
    generate_random_qip_model
using QIPresolve.ModelIO: save_moi

const DEFAULT_TARGET = joinpath("instances", "main_benchmark")

const RANDOM_COUNT = 30
const RANDOM_TYPES = ("bilinear", "separable", "pure", "generic")
const RANDOM_NVARS = 50
const RANDOM_NCONS = 50
const RANDOM_SEED_BASE = 10_000
const RANDOM_P_CON_EQ = 0.5
const RANDOM_VAR_THRESHOLD_LB = -10
const RANDOM_VAR_THRESHOLD_UB = 10
const RANDOM_P_VAR_IS_CANDIDATE = 0.2
const RANDOM_TERM_PROBABILITY = 0.2
const RANDOM_COEFF_LB = -50
const RANDOM_COEFF_UB = 50
const RANDOM_FORCE_DIAG_EVEN = false
const RANDOM_FORCE_LIN_EVEN = false
const RANDOM_FORCE_BILIN_EVEN = true
const RANDOM_FORCE_FEASIBILITY = true
const RANDOM_CONSTRAINT_SLACK_RANGE = collect(-5:5)
const RANDOM_CONSTRAINT_SLACK_RANGE_LABEL = "-5:5"

const EMBEDDING_COUNT = 15
const EMBEDDING_EXACTNESSES = ("exact", "inexact")
const EMBEDDING_GRAPH_TYPES = ("2con", "laman")
const EMBEDDING_ANCHORINGS = ("unanchored", "anchored")
const EMBEDDING_R = 100
const EMBEDDING_EDGE_DENSITY = 0.05
const EMBEDDING_NUM_ANCHORS = 3
const EMBEDDING_PH2 = 0.5
const EMBEDDING_MAX_COORD_TRIES = 10_000
const EMBEDDING_MAX_GLOBAL_TRIES = 10_000
const EMBEDDING_MAX_TRIES_H1 = 200
const EMBEDDING_MAX_TRIES_H2 = 300
const INFEASIBLE_EMBEDDING_COUNT = 15
const INFEASIBLE_EMBEDDING_STRATEGIES = (:bounding_box, :vertex_contraction)
const INFEASIBLE_EMBEDDING_SEED_BASE = 40_000
const INFEASIBLE_EMBEDDING_BOX_MARGIN = 1
const MAX_BOUNDING_BOX_CANDIDATE_TRIES = 10_000

function usage()
    return """
    Usage:
      julia --project=. scripts/generate_main_benchmark_instances.jl [--target PATH] [--force]

    Options:
      --target PATH   Output directory, default $DEFAULT_TARGET
      --force         Atomically replace an existing nonempty target
      -h, --help      Show this help
    """
end

function parse_args(args::Vector{String})
    target = DEFAULT_TARGET
    force = false
    index = 1
    while index <= length(args)
        arg = args[index]
        if arg in ("-h", "--help")
            println(usage())
            return nothing
        elseif startswith(arg, "--target=")
            target = split(arg, "="; limit = 2)[2]
        elseif arg == "--target"
            index += 1
            index <= length(args) || error("Missing value for --target")
            target = args[index]
        elseif arg == "--force"
            force = true
        else
            error("Unknown argument: $arg")
        end
        index += 1
    end
    isempty(strip(target)) && error("--target must be nonempty")
    return (target = abspath(target), force = force)
end

function random_term_probabilities(type::AbstractString)
    type == "bilinear" && return (bilin = RANDOM_TERM_PROBABILITY, diag = 0.0, lin = 0.0)
    type == "separable" && return (bilin = 0.0, diag = RANDOM_TERM_PROBABILITY, lin = RANDOM_TERM_PROBABILITY)
    type == "pure" && return (bilin = RANDOM_TERM_PROBABILITY, diag = RANDOM_TERM_PROBABILITY, lin = 0.0)
    type == "generic" && return (
        bilin = RANDOM_TERM_PROBABILITY,
        diag = RANDOM_TERM_PROBABILITY,
        lin = RANDOM_TERM_PROBABILITY,
    )
    error("Unknown random instance type: $type")
end

function random_specs(count::Int = RANDOM_COUNT)
    count >= 1 || error("random instance count must be positive")
    return [
        (
            type = type,
            replicate = replicate,
            seed = RANDOM_SEED_BASE + replicate - 1,
            probabilities = random_term_probabilities(type),
        )
        for type in RANDOM_TYPES for replicate in 1:count
    ]
end

function embedding_seed_base(exactness::AbstractString, graph_type::AbstractString)
    exactness == "exact" && graph_type == "2con" && return 20_000
    exactness == "exact" && graph_type == "laman" && return 21_000
    exactness == "inexact" && graph_type == "2con" && return 30_000
    exactness == "inexact" && graph_type == "laman" && return 31_000
    error("Unknown embedding family: $exactness/$graph_type")
end

function embedding_specs(count::Int = EMBEDDING_COUNT)
    count >= 1 || error("embedding instance count must be positive")
    feasible = [
        (
            exactness = exactness,
            graph_type = graph_type,
            anchoring = anchoring,
            feasibility = "feasible",
            infeas_strategy = nothing,
            infeas_base = nothing,
            replicate = replicate,
            seed = embedding_seed_base(exactness, graph_type) + replicate - 1,
            n = exactness == "exact" ? 26 : 28,
            alpha = exactness == "exact" ? 0.0 : 0.01,
            num_anchors = anchoring == "anchored" ? EMBEDDING_NUM_ANCHORS : 0,
        )
        for exactness in EMBEDDING_EXACTNESSES
        for graph_type in EMBEDDING_GRAPH_TYPES
        for anchoring in EMBEDDING_ANCHORINGS
        for replicate in 1:count
    ]

    infeasible = [
        (
            exactness = "exact",
            graph_type = "globally_rigid",
            anchoring = "anchored",
            feasibility = "infeasible",
            infeas_strategy = strategy,
            infeas_base = :globally_rigid,
            replicate = replicate,
            seed = INFEASIBLE_EMBEDDING_SEED_BASE + replicate - 1,
            n = 26,
            alpha = 0.0,
            num_anchors = EMBEDDING_NUM_ANCHORS,
        )
        for strategy in INFEASIBLE_EMBEDDING_STRATEGIES
        for replicate in 1:INFEASIBLE_EMBEDDING_COUNT
    ]

    return vcat(feasible, infeasible)
end

random_instance_name(spec) = "random_$(spec.type)_$(lpad(spec.replicate, 3, '0'))"

function embedding_instance_name(spec)
    suffix = lpad(spec.replicate, 3, '0')
    spec.feasibility == "infeasible" &&
        return "embedding_exact_infeasible_$(spec.infeas_strategy)_$suffix"
    return "embedding_$(spec.exactness)_$(spec.anchoring)_$(spec.graph_type)_$suffix"
end

function random_constraint_counts(model::JuMP.Model)
    scalar_functions = (JuMP.AffExpr, JuMP.QuadExpr)
    equalities = sum(
        length(JuMP.all_constraints(model, function_type, MOI.EqualTo{Float64}))
        for function_type in scalar_functions
    )
    inequalities = sum(
        length(JuMP.all_constraints(model, function_type, MOI.Interval{Float64}))
        for function_type in scalar_functions
    )
    return equalities, inequalities
end

function structural_constraint_count(model::JuMP.Model)
    return sum(
        length(JuMP.all_constraints(model, function_type, set_type))
        for (function_type, set_type) in JuMP.list_of_constraint_types(model)
        if function_type != JuMP.VariableRef
    )
end

function generate_random_instance!(directory::AbstractString, spec)
    probabilities = spec.probabilities
    model, _ = generate_random_qip_model(
        RANDOM_NVARS,
        RANDOM_NCONS;
        p_con_eq = RANDOM_P_CON_EQ,
        var_threshold_lb = RANDOM_VAR_THRESHOLD_LB,
        var_threshold_ub = RANDOM_VAR_THRESHOLD_UB,
        p_var_is_candidate = RANDOM_P_VAR_IS_CANDIDATE,
        p_var_bilin = probabilities.bilin,
        p_var_diag = probabilities.diag,
        p_var_lin = probabilities.lin,
        coeff_lb = RANDOM_COEFF_LB,
        coeff_ub = RANDOM_COEFF_UB,
        force_diag_even = RANDOM_FORCE_DIAG_EVEN,
        force_lin_even = RANDOM_FORCE_LIN_EVEN,
        force_bilin_even = RANDOM_FORCE_BILIN_EVEN,
        force_feasibility = RANDOM_FORCE_FEASIBILITY,
        constraint_slack_range = RANDOM_CONSTRAINT_SLACK_RANGE,
        seed = spec.seed,
    )

    instance_name = random_instance_name(spec)
    file_name = "$instance_name.lp"
    equalities, inequalities = random_constraint_counts(model)
    save_moi(JuMP.backend(model), joinpath(directory, file_name))

    return (
        instance_name = instance_name,
        file_name = file_name,
        subtype = spec.type,
        replicate = string(spec.replicate),
        seed = string(spec.seed),
        nvars = string(RANDOM_NVARS),
        ncons_requested = string(RANDOM_NCONS),
        ncons_realized = string(equalities + inequalities),
        eq_constraints = string(equalities),
        ineq_constraints = string(inequalities),
        p_con_eq = string(RANDOM_P_CON_EQ),
        var_threshold_lb = string(RANDOM_VAR_THRESHOLD_LB),
        var_threshold_ub = string(RANDOM_VAR_THRESHOLD_UB),
        p_var_is_candidate = string(RANDOM_P_VAR_IS_CANDIDATE),
        p_var_bilin = string(probabilities.bilin),
        p_var_diag = string(probabilities.diag),
        p_var_lin = string(probabilities.lin),
        coeff_lb = string(RANDOM_COEFF_LB),
        coeff_ub = string(RANDOM_COEFF_UB),
        force_diag_even = string(RANDOM_FORCE_DIAG_EVEN),
        force_lin_even = string(RANDOM_FORCE_LIN_EVEN),
        force_bilin_even = string(RANDOM_FORCE_BILIN_EVEN),
        force_feasibility = string(RANDOM_FORCE_FEASIBILITY),
        constraint_slack_range = RANDOM_CONSTRAINT_SLACK_RANGE_LABEL,
    )
end

function embedding_edge_count(spec)
    spec.graph_type == "laman" && return 2 * spec.n - 3
    spec.graph_type == "globally_rigid" && return 2 * spec.n - 2
    spec.graph_type == "2con" || error("Unknown embedding graph type: $(spec.graph_type)")
    minimum_edges = spec.n
    maximum_edges = spec.n * (spec.n - 1) ÷ 2
    return round(Int, minimum_edges + EMBEDDING_EDGE_DENSITY * (maximum_edges - minimum_edges))
end

function validate_model_variable_domains(model::JuMP.Model)
    conflicts = [
        (name = JuMP.name(var), lower = JuMP.lower_bound(var), upper = JuMP.upper_bound(var))
        for var in JuMP.all_variables(model)
        if JuMP.lower_bound(var) > JuMP.upper_bound(var)
    ]
    isempty(conflicts) || error("model contains contradictory variable bounds: $conflicts")
    return nothing
end

function generate_embedding_instance!(directory::AbstractString, spec; seed::Int = spec.seed)
    model, _, _ = if spec.feasibility == "infeasible"
        generate_likely_infeasible_embedding_instance(
            spec.n;
            strategy = spec.infeas_strategy,
            base = spec.infeas_base,
            seed = seed,
            num_anchors = spec.num_anchors,
            alpha = spec.alpha,
            box_margin = INFEASIBLE_EMBEDDING_BOX_MARGIN,
            contraction_vertices = nothing,
            R = EMBEDDING_R,
            max_global_tries = EMBEDDING_MAX_GLOBAL_TRIES,
            max_tries_H2 = EMBEDDING_MAX_TRIES_H2,
        )
    elseif spec.graph_type == "2con"
        generate_2_connected_instance(
            spec.n;
            R = EMBEDDING_R,
            edge_density = EMBEDDING_EDGE_DENSITY,
            seed = seed,
            max_coord_tries = EMBEDDING_MAX_COORD_TRIES,
            num_anchors = spec.num_anchors,
            alpha = spec.alpha,
        )
    elseif spec.graph_type == "laman"
        generate_laman_instance(
            spec.n;
            R = EMBEDDING_R,
            pH2 = EMBEDDING_PH2,
            seed = seed,
            max_global_tries = EMBEDDING_MAX_GLOBAL_TRIES,
            max_tries_H1 = EMBEDDING_MAX_TRIES_H1,
            max_tries_H2 = EMBEDDING_MAX_TRIES_H2,
            num_anchors = spec.num_anchors,
            alpha = spec.alpha,
        )
    else
        error("Unknown embedding graph type: $(spec.graph_type)")
    end

    validate_model_variable_domains(model)
    instance_name = embedding_instance_name(spec)
    file_name = "$instance_name.lp"
    save_moi(JuMP.backend(model), joinpath(directory, file_name))

    is_2con = spec.graph_type == "2con"
    is_laman = spec.graph_type == "laman"
    is_infeasible = spec.feasibility == "infeasible"
    return (
        instance_name = instance_name,
        file_name = file_name,
        exactness = spec.exactness,
        graph_type = spec.graph_type,
        anchoring = spec.anchoring,
        feasibility = spec.feasibility,
        replicate = string(spec.replicate),
        seed = string(seed),
        n = string(spec.n),
        nvars = string(JuMP.num_variables(model)),
        ncons_realized = string(structural_constraint_count(model)),
        edge_count = string(embedding_edge_count(spec)),
        R = string(EMBEDDING_R),
        num_anchors = string(spec.num_anchors),
        alpha = string(spec.alpha),
        edge_density = is_2con ? string(EMBEDDING_EDGE_DENSITY) : "",
        pH2 = is_laman ? string(EMBEDDING_PH2) : "",
        max_coord_tries = is_2con ? string(EMBEDDING_MAX_COORD_TRIES) : "",
        max_global_tries = is_2con ? "" : string(EMBEDDING_MAX_GLOBAL_TRIES),
        max_tries_H1 = is_laman ? string(EMBEDDING_MAX_TRIES_H1) : "",
        max_tries_H2 = is_2con ? "" : string(EMBEDDING_MAX_TRIES_H2),
        infeas_strategy = is_infeasible ? string(spec.infeas_strategy) : "",
        infeas_base = is_infeasible ? string(spec.infeas_base) : "",
        box_margin = spec.infeas_strategy == :bounding_box ? string(INFEASIBLE_EMBEDDING_BOX_MARGIN) : "",
        contraction_vertices = spec.infeas_strategy == :vertex_contraction ? "auto" : "",
    )
end


function generate_bounding_box_instance_with_retry!(
        directory::AbstractString, spec, initial_seed::Int
    )
    candidate_seed = initial_seed
    for _ in 1:MAX_BOUNDING_BOX_CANDIDATE_TRIES
        try
            return generate_embedding_instance!(directory, spec; seed = candidate_seed)
        catch err
            err isa BoundingBoxCandidateRejected || rethrow()
            println("rejected bounding-box seed $candidate_seed: $(sprint(showerror, err))")
            candidate_seed == typemax(Int) && error("cannot advance bounding-box candidate seed")
            candidate_seed += 1
        end
    end
    error(
        "failed to find a valid bounding-box candidate after " *
        "$MAX_BOUNDING_BOX_CANDIDATE_TRIES seeds starting at $initial_seed"
    )
end

function validate_specs(random, embedding)
    names = [random_instance_name(spec) for spec in random]
    append!(names, [embedding_instance_name(spec) for spec in embedding])
    allunique(names) || error("Benchmark instance names are not unique")
    return nothing
end

function validate_staging_directory(directory::AbstractString, random_rows, embedding_rows)
    expected_files = Set(row.file_name for row in Iterators.flatten((random_rows, embedding_rows)))
    actual_files = Set(filter(endswith(".lp"), readdir(directory)))
    actual_files == expected_files || error("Generated LP files do not match the benchmark specification")

    random_csv_rows = collect(CSV.File(joinpath(directory, "random_instances.csv")))
    embedding_csv_rows = collect(CSV.File(joinpath(directory, "embedding_instances.csv")))
    length(random_csv_rows) == length(random_rows) || error("Random CSV row count is incorrect")
    length(embedding_csv_rows) == length(embedding_rows) || error("Embedding CSV row count is incorrect")
    return nothing
end

function ensure_target_available(target::AbstractString; force::Bool = false)
    (ispath(target) || islink(target)) || return nothing
    isdir(target) && !islink(target) || error("Target exists and is not a directory: $target")
    isempty(readdir(target)) || force || error(
        "Refusing to overwrite nonempty target directory without --force: $target"
    )
    return nothing
end

function publish_staging_directory!(
        staging::AbstractString,
        target::AbstractString;
        force::Bool = false,
    )
    ensure_target_available(target; force = force)
    if !isdir(target)
        mv(staging, target)
        return nothing
    elseif isempty(readdir(target))
        rm(target)
        mv(staging, target)
        return nothing
    end

    backup = "$target.backup"
    (ispath(backup) || islink(backup)) && error("Backup path already exists: $backup")
    mv(target, backup)
    try
        mv(staging, target)
    catch
        (ispath(target) || islink(target)) || mv(backup, target)
        rethrow()
    end
    rm(backup; recursive = true)
    return nothing
end

function generate_dataset(
        target::AbstractString;
        random = random_specs(),
        embedding = embedding_specs(),
        force::Bool = false,
    )
    target = abspath(target)
    ensure_target_available(target; force = force)
    validate_specs(random, embedding)
    mkpath(dirname(target))
    staging = mktempdir(dirname(target); prefix = "$(basename(target)).staging-")

    try
        random_rows = NamedTuple[]
        for spec in random
            row = generate_random_instance!(staging, spec)
            push!(random_rows, row)
            println("saved $(row.instance_name)")
        end

        embedding_rows = NamedTuple[]
        next_bounding_box_seed = nothing
        for spec in embedding
            row = if spec.infeas_strategy == :bounding_box
                initial_seed = next_bounding_box_seed === nothing ?
                    spec.seed : max(spec.seed, next_bounding_box_seed)
                bounding_row = generate_bounding_box_instance_with_retry!(
                    staging, spec, initial_seed
                )
                next_bounding_box_seed = parse(Int, bounding_row.seed) + 1
                bounding_row
            else
                generate_embedding_instance!(staging, spec)
            end
            push!(embedding_rows, row)
            println("saved $(row.instance_name)")
        end

        CSV.write(joinpath(staging, "random_instances.csv"), random_rows)
        CSV.write(joinpath(staging, "embedding_instances.csv"), embedding_rows)
        validate_staging_directory(staging, random_rows, embedding_rows)
        publish_staging_directory!(staging, target; force = force)
        staging = ""
    finally
        !isempty(staging) && isdir(staging) && rm(staging; recursive = true)
    end

    println("generated $(length(random)) random and $(length(embedding)) embedding instances in $target")
    return target
end

function main(args::Vector{String})
    config = parse_args(args)
    config === nothing && return nothing
    return generate_dataset(config.target; force = config.force)
end

end # module MainBenchmarkInstanceGenerator

if _RUNNING_AS_SCRIPT
    MainBenchmarkInstanceGenerator.main(ARGS)
end
