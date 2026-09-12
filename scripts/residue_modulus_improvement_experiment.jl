#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module ResidueModulusImprovementExperiment

using CSV
using JuMP: backend
using Printf

import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC
using QIPresolve.InstanceGeneration: generate_random_qip_model

const DEFAULT_COUNT = 20
const DEFAULT_N_VARS = 50
const DEFAULT_N_CONSTRAINTS = 50
const DEFAULT_SEED_BASE = 20000
const DEFAULT_SEED_STEP = 1
const DEFAULT_MAX_MODULUS = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD
const DEFAULT_CONSTRAINT_SLACK_RANGE = [collect(-100:-1); collect(1:100)]
const DEFAULT_CONSTRAINT_SLACK_RANGE_LABEL = "-100:-1,1:100"

const CLI_KEYS = Dict(
    "count" => :count,
    "nvars" => :nvars,
    "n-vars" => :nvars,
    "ncons" => :ncons,
    "n-cons" => :ncons,
    "nconstraints" => :ncons,
    "n-constraints" => :ncons,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "max-modulus" => :max_modulus,
    "max-moduli" => :max_modulus,
    "m" => :max_modulus,
    "residue-threshold" => :max_modulus,
    "treewidth-threshold" => :treewidth_threshold,
    "output" => :output_path,
    "p-con-eq" => :p_con_eq,
    "var-threshold-lb" => :var_threshold_lb,
    "var-threshold-ub" => :var_threshold_ub,
    "p-var-is-candidate" => :p_var_is_candidate,
    "p-var-candidate" => :p_var_is_candidate,
    "p-var-bilin" => :p_var_bilin,
    "p-var-diag" => :p_var_diag,
    "p-var-lin" => :p_var_lin,
    "coeff-lb" => :coeff_lb,
    "coeff-ub" => :coeff_ub,
    "force-diag-even" => :force_diag_even,
    "force-lin-even" => :force_lin_even,
    "force-feasibility" => :force_feasibility,
    "constraint-slack-range" => :constraint_slack_range,
)

Base.@kwdef struct CliConfig
    count::Int = DEFAULT_COUNT
    nvars::Int = DEFAULT_N_VARS
    ncons::Int = DEFAULT_N_CONSTRAINTS
    seed_base::Int = DEFAULT_SEED_BASE
    seed_step::Int = DEFAULT_SEED_STEP
    max_modulus::Int = DEFAULT_MAX_MODULUS
    treewidth_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD
    output_path::Union{Nothing, String} = nothing
    p_con_eq::Float64 = 0.0
    var_threshold_lb::Int = -10
    var_threshold_ub::Int = 10
    p_var_is_candidate::Float64 = 0.1
    p_var_bilin::Float64 = 0.3
    p_var_diag::Float64 = 0.3
    p_var_lin::Float64 = 0.3
    coeff_lb::Int = -20
    coeff_ub::Int = 20
    force_diag_even::Bool = false
    force_lin_even::Bool = false
    force_feasibility::Bool = true
    constraint_slack_range::Vector{Int} = copy(DEFAULT_CONSTRAINT_SLACK_RANGE)
    constraint_slack_range_label::String = DEFAULT_CONSTRAINT_SLACK_RANGE_LABEL
end

struct StrategySpec
    name::String
    moduli::Vector{Int}
end

Base.@kwdef mutable struct ModulusStepCounts
    modulus::Int
    modulus_order::Int
    evaluated_constraints::Int = 0
    bound_tightened_constraints::Int = 0
    lower_bound_tightenings::Int = 0
    upper_bound_tightenings::Int = 0
    infeasibilities::Int = 0
end

Base.@kwdef mutable struct StrategyResult
    name::String
    moduli::Vector{Int}
    modulus_counts::Vector{ModulusStepCounts}
    processed_constraints::Int = 0
    improved_constraints::Int = 0
    infeasibilities::Int = 0
    total_relative_bound_tightening::Float64 = 0.0
    strategy_time_sec::Float64 = 0.0
end

function StrategyResult(spec::StrategySpec)
    counts = [
        ModulusStepCounts(modulus = modulus, modulus_order = index)
        for (index, modulus) in enumerate(spec.moduli)
    ]
    return StrategyResult(name = spec.name, moduli = spec.moduli, modulus_counts = counts)
end

function usage()
    return """
    Usage:
      julia --project=. scripts/residue_modulus_improvement_experiment.jl [options]

    Options:
      --count n                       Random instances to generate, default $DEFAULT_COUNT
      --nvars n                       Variables per instance, default $DEFAULT_N_VARS
      --ncons n                       Constraints per instance, default $DEFAULT_N_CONSTRAINTS
      --seed-base n                   First random seed, default $DEFAULT_SEED_BASE
      --seed-step n                   Seed increment, default $DEFAULT_SEED_STEP
      --max-modulus n                 Maximum residue modulus M, default $DEFAULT_MAX_MODULUS
      --residue-threshold n           Backward-compatible alias for --max-modulus
      --treewidth-threshold n         Residue DP treewidth threshold
      --output path                   Optional CSV output path
      --p-var-is-candidate p          Candidate variable probability, default 0.1
      --p-var-bilin p                 Bilinear term probability among candidates, default 0.3
      --p-var-diag p                  Diagonal term probability among candidates, default 0.3
      --p-var-lin p                   Linear term probability among candidates, default 0.3
      --constraint-slack-range spec   Comma-separated ints/ranges, default $DEFAULT_CONSTRAINT_SLACK_RANGE_LABEL
      -h, --help                      Show this help
    """
end

function parse_int(value::AbstractString, name::AbstractString)::Int
    parsed = tryparse(Int, value)
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function parse_float(value::AbstractString, name::AbstractString)::Float64
    parsed = tryparse(Float64, value)
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function parse_bool(value::AbstractString, name::AbstractString)::Bool
    normalized = lowercase(strip(value))
    normalized in ("true", "1", "yes") && return true
    normalized in ("false", "0", "no") && return false
    error("Invalid $name: $value. Expected true/false, 1/0, or yes/no.")
end

function parse_constraint_slack_range(value::AbstractString)
    normalized = replace(strip(value), " " => "")
    isempty(normalized) && error("Invalid constraint slack range: value must be nonempty")

    values = Int[]
    for part in split(normalized, ","; keepempty = false)
        if occursin(":", part)
            bounds = split(part, ":"; limit = 2)
            length(bounds) == 2 || error("Invalid constraint slack range part: $part")
            lower = parse_int(bounds[1], "constraint slack lower bound")
            upper = parse_int(bounds[2], "constraint slack upper bound")
            lower <= upper || error("constraint slack range lower bound must be <= upper bound")
            append!(values, lower:upper)
        else
            push!(values, parse_int(part, "constraint slack value"))
        end
    end

    isempty(values) && error("Invalid constraint slack range: value must be nonempty")
    return values, normalized
end

function _option_value(args::Vector{String}, index::Int, option::String)
    if occursin("=", option)
        key, value = split(option, "="; limit = 2)
        return key, value, 1
    end

    index < length(args) || error("Missing value for option $option")
    return option, args[index + 1], 2
end

function parse_raw_options(args::Vector{String})
    options = Dict{Symbol, String}()
    index = 1
    while index <= length(args)
        arg = args[index]
        if arg in ("-h", "--help")
            println(usage())
            return nothing
        end

        startswith(arg, "--") || error("Unexpected positional argument: $arg")
        raw_key, value, consumed = _option_value(args, index, arg)
        lookup_key = lowercase(replace(raw_key[3:end], "_" => "-"))
        key = get(CLI_KEYS, lookup_key, nothing)
        key === nothing && error("Unknown option: $raw_key")
        options[key] = value
        index += consumed
    end
    return options
end

function validate_probability(name::AbstractString, value::Float64)
    isfinite(value) && 0.0 <= value <= 1.0 ||
        error("$name must be in [0, 1], got $value")
    return value
end

function build_config(args::Vector{String})::Union{Nothing, CliConfig}
    options = parse_raw_options(args)
    options === nothing && return nothing

    slack_range, slack_label = haskey(options, :constraint_slack_range) ?
        parse_constraint_slack_range(options[:constraint_slack_range]) :
        (copy(DEFAULT_CONSTRAINT_SLACK_RANGE), DEFAULT_CONSTRAINT_SLACK_RANGE_LABEL)

    config = CliConfig(
        count = parse_int(get(options, :count, string(DEFAULT_COUNT)), "count"),
        nvars = parse_int(get(options, :nvars, string(DEFAULT_N_VARS)), "nvars"),
        ncons = parse_int(get(options, :ncons, string(DEFAULT_N_CONSTRAINTS)), "ncons"),
        seed_base = parse_int(get(options, :seed_base, string(DEFAULT_SEED_BASE)), "seed_base"),
        seed_step = parse_int(get(options, :seed_step, string(DEFAULT_SEED_STEP)), "seed_step"),
        max_modulus = parse_int(
            get(options, :max_modulus, string(DEFAULT_MAX_MODULUS)),
            "max_modulus",
        ),
        treewidth_threshold = parse_int(
            get(options, :treewidth_threshold, string(QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD)),
            "treewidth_threshold",
        ),
        output_path = haskey(options, :output_path) ? abspath(options[:output_path]) : nothing,
        p_con_eq = parse_float(get(options, :p_con_eq, "0.0"), "p_con_eq"),
        var_threshold_lb = parse_int(get(options, :var_threshold_lb, "-10"), "var_threshold_lb"),
        var_threshold_ub = parse_int(get(options, :var_threshold_ub, "10"), "var_threshold_ub"),
        p_var_is_candidate = parse_float(get(options, :p_var_is_candidate, "0.1"), "p_var_is_candidate"),
        p_var_bilin = parse_float(get(options, :p_var_bilin, "0.3"), "p_var_bilin"),
        p_var_diag = parse_float(get(options, :p_var_diag, "0.3"), "p_var_diag"),
        p_var_lin = parse_float(get(options, :p_var_lin, "0.3"), "p_var_lin"),
        coeff_lb = parse_int(get(options, :coeff_lb, "-50"), "coeff_lb"),
        coeff_ub = parse_int(get(options, :coeff_ub, "50"), "coeff_ub"),
        force_diag_even = parse_bool(get(options, :force_diag_even, "false"), "force_diag_even"),
        force_lin_even = parse_bool(get(options, :force_lin_even, "false"), "force_lin_even"),
        force_feasibility = parse_bool(get(options, :force_feasibility, "true"), "force_feasibility"),
        constraint_slack_range = slack_range,
        constraint_slack_range_label = slack_label,
    )

    config.count >= 1 || error("count must be >= 1")
    config.nvars >= 1 || error("nvars must be >= 1")
    config.ncons >= 0 || error("ncons must be >= 0")
    config.seed_step >= 0 || error("seed_step must be >= 0")
    config.max_modulus >= 2 || error("max_modulus must be >= 2")
    config.treewidth_threshold >= 0 || error("treewidth_threshold must be >= 0")
    config.var_threshold_lb <= config.var_threshold_ub ||
        error("var_threshold_lb must be <= var_threshold_ub")
    config.coeff_lb <= config.coeff_ub || error("coeff_lb must be <= coeff_ub")

    validate_probability("p_con_eq", config.p_con_eq)
    validate_probability("p_var_is_candidate", config.p_var_is_candidate)
    validate_probability("p_var_bilin", config.p_var_bilin)
    validate_probability("p_var_diag", config.p_var_diag)
    validate_probability("p_var_lin", config.p_var_lin)

    return config
end

function random_qip_kwargs(config::CliConfig)
    return (
        p_con_eq = config.p_con_eq,
        var_threshold_lb = config.var_threshold_lb,
        var_threshold_ub = config.var_threshold_ub,
        p_var_is_candidate = config.p_var_is_candidate,
        p_var_bilin = config.p_var_bilin,
        p_var_diag = config.p_var_diag,
        p_var_lin = config.p_var_lin,
        coeff_lb = config.coeff_lb,
        coeff_ub = config.coeff_ub,
        force_diag_even = config.force_diag_even,
        force_lin_even = config.force_lin_even,
        force_feasibility = config.force_feasibility,
        constraint_slack_range = config.constraint_slack_range,
    )
end

function build_random_qip_core_model(config::CliConfig, seed::Int)
    jump_model, _ = generate_random_qip_model(
        config.nvars,
        config.ncons;
        random_qip_kwargs(config)...,
        seed = seed,
    )
    return QIP.build_model(QIP.from_moi(backend(jump_model)))
end

function one_constraint_model(model::PC.QPModel, con::PC.Constraint)
    return PC.QPModel(model.vars, [con], model.obj_expr, model.obj_sense)
end

function is_prime(candidate::Int)
    candidate < 2 && return false
    candidate == 2 && return true
    iseven(candidate) && return false

    divisor = 3
    while divisor <= div(candidate, divisor)
        candidate % divisor == 0 && return false
        divisor += 2
    end
    return true
end

function prime_moduli(upper::Int)
    upper < 2 && return Int[]
    return [candidate for candidate in 2:upper if is_prime(candidate)]
end

function prime_power_moduli(upper::Int)
    upper < 2 && return Int[]

    values = Set{Int}()
    for prime in prime_moduli(upper)
        value = prime
        while value <= upper
            push!(values, value)
            value > div(upper, prime) && break
            value *= prime
        end
    end

    return sort!(collect(values))
end

function strategy_specs(config::CliConfig)
    upper = config.max_modulus

    return [
        StrategySpec("all_2_M", collect(upper:-1:cld(upper, 2))),
        StrategySpec("primes_2_M", prime_moduli(upper)),
        StrategySpec("prime_powers_2_M", prime_power_moduli(upper)),
    ]
end

bound_snapshot(con::PC.Constraint) = (
    lhs = con.lhs,
    rhs = con.rhs,
    equality = PC.is_equality(con),
)

function relative_bound_tightening(before, con::PC.Constraint)
    isfinite(before.lhs) && isfinite(before.rhs) || return 0.0

    baseline_range = before.rhs - before.lhs
    baseline_range > 0.0 || return 0.0

    improvement = 0.0
    isfinite(con.lhs) && (improvement += max(0.0, con.lhs - before.lhs))
    isfinite(con.rhs) && (improvement += max(0.0, before.rhs - con.rhs))

    return improvement / baseline_range
end

function run_sequential_strategy!(
        model::PC.QPModel,
        con::PC.Constraint,
        result::StrategyResult,
        treewidth_threshold::Int,
    )
    model.infeasible && return true

    status = PC._residue_constraint_status(con, model.vars)
    if status == :infeasible
        model.infeasible = true
        return true
    elseif status == :skip || isempty(result.moduli)
        return false
    end

    standardized = PC._standardize_residue_constraint(con, model.vars, model.obj_expr)
    cache = PC._ResidueCacheEntry[]

    for counts in result.modulus_counts
        modulus = counts.modulus
        counts.evaluated_constraints += 1
        residue_result = PC._compute_achievable_residues(
            modulus,
            standardized.con,
            standardized.var_bounds;
            treewidth_threshold = treewidth_threshold,
        )

        residue_result.saturated && continue

        push!(cache, PC._ResidueCacheEntry(modulus, residue_result.residues))
        before = bound_snapshot(con)
        PC._reapply_residue_cache!(con, cache, standardized.constraint_shift)

        lower_tightened = isfinite(before.lhs) && con.lhs > before.lhs
        upper_tightened = isfinite(before.rhs) && con.rhs < before.rhs
        if lower_tightened || upper_tightened
            counts.bound_tightened_constraints += 1
            lower_tightened && (counts.lower_bound_tightenings += 1)
            upper_tightened && (counts.upper_bound_tightenings += 1)
        end

        if con.lhs > con.rhs
            model.infeasible = true
            counts.infeasibilities += 1
            return true
        end
    end

    return false
end

function record_strategy_trial!(
        result::StrategyResult,
        model::PC.QPModel,
        con::PC.Constraint,
        treewidth_threshold::Int,
    )
    trial_con = deepcopy(con)
    trial_model = one_constraint_model(model, trial_con)
    before = bound_snapshot(trial_con)

    result.processed_constraints += 1
    start_time = time()
    infeasible = run_sequential_strategy!(
        trial_model,
        trial_con,
        result,
        treewidth_threshold,
    )
    result.strategy_time_sec += time() - start_time
    relative_tightening = relative_bound_tightening(before, trial_con)

    result.total_relative_bound_tightening += relative_tightening
    relative_tightening > 0.0 && (result.improved_constraints += 1)
    (infeasible || trial_con.lhs > trial_con.rhs) && (result.infeasibilities += 1)

    return result
end

function record_strategy_trials!(
        results::Vector{StrategyResult},
        model::PC.QPModel,
        con::PC.Constraint,
        treewidth_threshold::Int,
    )
    for result in results
        record_strategy_trial!(result, model, con, treewidth_threshold)
    end
    return results
end

rate(numerator::Real, denominator::Int) = denominator == 0 ? 0.0 : numerator / denominator

function row(result::StrategyResult)
    return (
        strategy = result.name,
        moduli = join(result.moduli, " "),
        num_moduli = length(result.moduli),
        processed_constraints = result.processed_constraints,
        improved_constraints = result.improved_constraints,
        avg_relative_bound_tightening = rate(
            result.total_relative_bound_tightening,
            result.processed_constraints,
        ),
        strategy_time_sec = result.strategy_time_sec,
        infeasibilities = result.infeasibilities,
    )
end

function detailed_row(config::CliConfig, result::StrategyResult, counts::ModulusStepCounts)
    summary = row(result)
    return (
        strategy = result.name,
        max_modulus = config.max_modulus,
        modulus_order = counts.modulus_order,
        modulus = counts.modulus,
        num_moduli = summary.num_moduli,
        processed_constraints = summary.processed_constraints,
        improved_constraints = summary.improved_constraints,
        avg_relative_bound_tightening = summary.avg_relative_bound_tightening,
        strategy_time_sec = summary.strategy_time_sec,
        infeasibilities = summary.infeasibilities,
        modulus_evaluated_constraints = counts.evaluated_constraints,
        modulus_bound_tightened_constraints = counts.bound_tightened_constraints,
        modulus_lower_bound_tightenings = counts.lower_bound_tightenings,
        modulus_upper_bound_tightenings = counts.upper_bound_tightenings,
        modulus_infeasibilities = counts.infeasibilities,
    )
end

function detailed_rows(config::CliConfig, results::Vector{StrategyResult})
    rows = NamedTuple[]
    for result in results
        for counts in result.modulus_counts
            push!(rows, detailed_row(config, result, counts))
        end
    end
    return rows
end

function run_experiment(config::CliConfig)
    strategies = strategy_specs(config)
    results = StrategyResult.(strategies)
    generated_instances = 0
    skipped_infeasible_instances = 0
    processed_constraints = 0

    for offset in 0:(config.count - 1)
        seed = config.seed_base + offset * config.seed_step
        model = build_random_qip_core_model(config, seed)
        generated_instances += 1
        PC.normalize!(model, nothing)

        if model.infeasible
            skipped_infeasible_instances += 1
            continue
        end

        for con in model.cons
            processed_constraints += 1
            record_strategy_trials!(
                results,
                model,
                con,
                config.treewidth_threshold,
            )
        end
    end

    rows = [row(result) for result in results]
    return (
        config = config,
        strategies = strategies,
        rows = rows,
        csv_rows = detailed_rows(config, results),
        generated_instances = generated_instances,
        skipped_infeasible_instances = skipped_infeasible_instances,
        processed_constraints = processed_constraints,
    )
end

function write_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    CSV.write(path, rows)
    return path
end

function print_config(result)
    config = result.config
    println("Residue modulus improvement experiment")
    println("instances = $(config.count)")
    println("nvars = $(config.nvars)")
    println("nconstraints = $(config.ncons)")
    println("seed_base = $(config.seed_base)")
    println("seed_step = $(config.seed_step)")
    println("max_modulus = $(config.max_modulus)")
    println("treewidth_threshold = $(config.treewidth_threshold)")
    println("p_var_is_candidate = $(config.p_var_is_candidate)")
    println("p_var_bilin = $(config.p_var_bilin)")
    println("p_var_diag = $(config.p_var_diag)")
    println("p_var_lin = $(config.p_var_lin)")
    println("constraint_slack_range = $(config.constraint_slack_range_label)")
    println("generated_instances = $(result.generated_instances)")
    println("skipped_infeasible_instances = $(result.skipped_infeasible_instances)")
    println("processed_constraints = $(result.processed_constraints)")
    println()
    return nothing
end

function print_table(rows)
    @printf(
        "%-20s %10s %20s %18s %30s %12s %15s\n",
        "strategy",
        "n_moduli",
        "processed_constraints",
        "improved_constraints",
        "avg_relative_bound_tightening",
        "time_sec",
        "infeasibilities",
    )
    for row in rows
        @printf(
            "%-20s %10d %20d %18d %30.6f %12.6f %15d\n",
            row.strategy,
            row.num_moduli,
            row.processed_constraints,
            row.improved_constraints,
            row.avg_relative_bound_tightening,
            row.strategy_time_sec,
            row.infeasibilities,
        )
    end
    return nothing
end

function main(args::Vector{String} = copy(ARGS))
    config = build_config(args)
    config === nothing && return nothing

    result = run_experiment(config)
    print_config(result)
    print_table(result.rows)

    if config.output_path !== nothing
        write_csv(config.output_path, result.csv_rows)
        println()
        println("Wrote CSV to $(config.output_path)")
    end

    return result
end

end # module

if _RUNNING_AS_SCRIPT
    ResidueModulusImprovementExperiment.main(copy(ARGS))
end
