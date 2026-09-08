#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module ResidueRandomQuadraticConstraintExperiment

using Printf
using Random

include("residue_constraint_experiment_common.jl")

const Common = ResidueConstraintExperimentCommon
const PC = Common.PC
const QuadTerm = Common.QuadTerm
const LinTerm = Common.LinTerm

const DEFAULT_COUNT = 1000
const DEFAULT_NVARS = (10,)
const DEFAULT_SEED_BASE = 31000
const DEFAULT_SEED_STEP = 1
const DEFAULT_DOMAIN_LB = 0
const DEFAULT_DOMAIN_UB = 5
const DEFAULT_DENSITY = 0.1
const DEFAULT_COEFF_LB = -50
const DEFAULT_COEFF_UB = 50
const DEFAULT_MAX_DISTINCT_COEFFS = 50
const DEFAULT_OFFSET_LB = 1
const DEFAULT_OFFSET_UB = 10
const DEFAULT_MAX_GENERATION_TRIES = 100
const DEFAULT_EXACT_ENUMERATION = false
const DEFAULT_MODULI_STRATEGY = "primes_lt_128"

const CLI_KEYS = Dict(
    "count" => :count,
    "nvars" => :nvars,
    "n-vars" => :nvars,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "density" => :density,
    "domain-lb" => :domain_lb,
    "domain-ub" => :domain_ub,
    "coeff-lb" => :coeff_lb,
    "coeff-ub" => :coeff_ub,
    "max-distinct-coeffs" => :max_distinct_coeffs,
    "max-distinct-coefficients" => :max_distinct_coeffs,
    "offset-lb" => :offset_lb,
    "offset-min" => :offset_lb,
    "bound-offset-lb" => :offset_lb,
    "offset-ub" => :offset_ub,
    "offset-max" => :offset_ub,
    "bound-offset-ub" => :offset_ub,
    "treewidth-threshold" => :treewidth_threshold,
    "max-generation-tries" => :max_generation_tries,
    "generation-tries" => :max_generation_tries,
    "exact-enumeration" => :exact_enumeration,
    "brute-force-enumeration" => :exact_enumeration,
    "bruteforce-enumeration" => :exact_enumeration,
    "moduli-strategy" => :moduli_strategies,
    "moduli-strategies" => :moduli_strategies,
    "strategy" => :moduli_strategies,
    "output" => :output_path,
)

Base.@kwdef struct CliConfig
    count::Int = DEFAULT_COUNT
    nvars::Vector{Int} = collect(DEFAULT_NVARS)
    seed_base::Int = DEFAULT_SEED_BASE
    seed_step::Int = DEFAULT_SEED_STEP
    density::Float64 = DEFAULT_DENSITY
    domain_lb::Int = DEFAULT_DOMAIN_LB
    domain_ub::Int = DEFAULT_DOMAIN_UB
    coeff_lb::Int = DEFAULT_COEFF_LB
    coeff_ub::Int = DEFAULT_COEFF_UB
    max_distinct_coeffs::Int = DEFAULT_MAX_DISTINCT_COEFFS
    offset_lb::Int = DEFAULT_OFFSET_LB
    offset_ub::Int = DEFAULT_OFFSET_UB
    treewidth_threshold::Int = Common.DEFAULT_TREEWIDTH_THRESHOLD
    max_generation_tries::Int = DEFAULT_MAX_GENERATION_TRIES
    exact_enumeration::Bool = DEFAULT_EXACT_ENUMERATION
    moduli_strategies::Vector{String} = String[DEFAULT_MODULI_STRATEGY]
    output_path::Union{Nothing, String} = nothing
end

struct ConstraintSample
    seed::Int
    nvars::Int
    model::PC.QPModel
    con::PC.Constraint
    x_star::Vector{Int}
    generation_tries::Int
end

Base.@kwdef mutable struct StrategyResult
    nvars::Int
    density::Float64
    domain_lb::Int
    domain_ub::Int
    max_distinct_coeffs::Int
    name::String
    moduli::Vector{Int}
    exact_enumeration::Bool = DEFAULT_EXACT_ENUMERATION
    constraints::Int = 0
    bounds_tightened::Int = 0
    constraints_tightened_to_equality::Int = 0
    total_relative_bound_range_reduction::Float64 = 0.0
    bounds_fully_tightened_to_optimal::Int = 0
    total_bound_gap_to_optimal::Float64 = 0.0
    total_residue_time_sec::Float64 = 0.0
    exact_assignments_per_constraint::Union{Missing, Int} = missing
end

function usage()
    return """
    Usage:
      julia --project=. scripts/residue_random_quadratic_constraint_experiment.jl [options]

    Options:
      --count n                       Constraints per n, default $DEFAULT_COUNT
      --nvars list                    Comma-separated n values, default $(join(DEFAULT_NVARS, ","))
      --seed-base n                   First random seed, default $DEFAULT_SEED_BASE
      --seed-step n                   Seed increment, default $DEFAULT_SEED_STEP
      --density p                     Probability for each coefficient to be nonzero, default $DEFAULT_DENSITY
      --domain-lb n                   Variable lower bound, default $DEFAULT_DOMAIN_LB
      --domain-ub n                   Variable upper bound, default $DEFAULT_DOMAIN_UB
      --coeff-lb n                    Coefficient lower bound, default $DEFAULT_COEFF_LB
      --coeff-ub n                    Coefficient upper bound, default $DEFAULT_COEFF_UB
      --max-distinct-coeffs n         Max distinct coefficients per generated constraint, default $DEFAULT_MAX_DISTINCT_COEFFS
      --offset-lb n                   Lower offset bound for sampled constraint slack, default $DEFAULT_OFFSET_LB
      --offset-ub n                   Upper offset bound for sampled constraint slack, default $DEFAULT_OFFSET_UB
      --treewidth-threshold n         Residue DP treewidth threshold
      --max-generation-tries n        Attempts to avoid degenerate generated constraints, default $DEFAULT_MAX_GENERATION_TRIES
      --exact-enumeration bool        Run brute-force exact bound comparison, default $DEFAULT_EXACT_ENUMERATION
      --moduli-strategy list          Comma-separated strategies: $(join(Common.moduli_family_strategy_names(), ", "))
      --output path                   Optional CSV output path
      -h, --help                      Show this help
    """
end

function _lookup_option_key(raw_key::AbstractString)
    startswith(raw_key, "--") || error("Unexpected positional argument: $raw_key")
    lookup_key = lowercase(replace(raw_key[3:end], "_" => "-"))
    key = get(CLI_KEYS, lookup_key, nothing)
    key === nothing && error("Unknown option: $raw_key")
    return key
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

        raw_key = arg
        value = ""
        consumed = 1

        if occursin("=", arg)
            raw_key, value = split(arg, "="; limit = 2)
            consumed = 1
        else
            _lookup_option_key(raw_key)
            index < length(args) || error("Missing value for option $raw_key")
            value = args[index + 1]
            consumed = 2
        end

        key = _lookup_option_key(raw_key)
        options[key] = value
        index += consumed
    end

    return options
end

coefficient_values(config::CliConfig) =
    Common.coefficient_values(config.coeff_lb, config.coeff_ub)

function coefficient_palette(rng::AbstractRNG, config::CliConfig)
    return Common.coefficient_palette(
        rng,
        config.coeff_lb,
        config.coeff_ub,
        config.max_distinct_coeffs,
    )
end

function validate_config(config::CliConfig)
    config.count >= 1 || error("count must be >= 1")
    isempty(config.nvars) && error("nvars must contain at least one value")
    all(>=(1), config.nvars) || error("all nvars values must be >= 1")
    config.seed_step >= 0 || error("seed_step must be >= 0")
    Common.validate_probability("density", config.density)
    config.domain_lb <= config.domain_ub || error("domain_lb must be <= domain_ub")
    config.coeff_lb <= config.coeff_ub || error("coeff_lb must be <= coeff_ub")
    coefficient_values(config)
    config.max_distinct_coeffs >= 1 || error("max_distinct_coeffs must be >= 1")
    config.offset_lb >= 0 || error("offset_lb must be >= 0")
    config.offset_lb <= config.offset_ub || error("offset_lb must be <= offset_ub")
    config.treewidth_threshold >= 0 || error("treewidth_threshold must be >= 0")
    config.max_generation_tries >= 1 || error("max_generation_tries must be >= 1")
    Common.normalize_moduli_strategies(config.moduli_strategies)

    return config
end

function build_config(args::Vector{String})::Union{Nothing, CliConfig}
    options = parse_raw_options(args)
    options === nothing && return nothing

    config = CliConfig(
        count = Common.parse_int(get(options, :count, string(DEFAULT_COUNT)), "count"),
        nvars = haskey(options, :nvars) ?
            Common.parse_int_list(options[:nvars], "nvars") :
            collect(DEFAULT_NVARS),
        seed_base = Common.parse_int(get(options, :seed_base, string(DEFAULT_SEED_BASE)), "seed_base"),
        seed_step = Common.parse_int(get(options, :seed_step, string(DEFAULT_SEED_STEP)), "seed_step"),
        density = Common.parse_float(
            get(options, :density, string(DEFAULT_DENSITY)),
            "density",
        ),
        domain_lb = Common.parse_int(get(options, :domain_lb, string(DEFAULT_DOMAIN_LB)), "domain_lb"),
        domain_ub = Common.parse_int(get(options, :domain_ub, string(DEFAULT_DOMAIN_UB)), "domain_ub"),
        coeff_lb = Common.parse_int(get(options, :coeff_lb, string(DEFAULT_COEFF_LB)), "coeff_lb"),
        coeff_ub = Common.parse_int(get(options, :coeff_ub, string(DEFAULT_COEFF_UB)), "coeff_ub"),
        max_distinct_coeffs = haskey(options, :max_distinct_coeffs) ?
            Common.parse_int(options[:max_distinct_coeffs], "max_distinct_coeffs") :
            DEFAULT_MAX_DISTINCT_COEFFS,
        offset_lb = Common.parse_int(get(options, :offset_lb, string(DEFAULT_OFFSET_LB)), "offset_lb"),
        offset_ub = Common.parse_int(get(options, :offset_ub, string(DEFAULT_OFFSET_UB)), "offset_ub"),
        treewidth_threshold = Common.parse_int(
            get(options, :treewidth_threshold, string(Common.DEFAULT_TREEWIDTH_THRESHOLD)),
            "treewidth_threshold",
        ),
        max_generation_tries = Common.parse_int(
            get(options, :max_generation_tries, string(DEFAULT_MAX_GENERATION_TRIES)),
            "max_generation_tries",
        ),
        exact_enumeration = Common.parse_bool(
            get(options, :exact_enumeration, string(DEFAULT_EXACT_ENUMERATION)),
            "exact_enumeration",
        ),
        moduli_strategies = haskey(options, :moduli_strategies) ?
            Common.normalize_moduli_strategies(
                Common.parse_string_list(options[:moduli_strategies], "moduli_strategies"),
            ) :
            String[DEFAULT_MODULI_STRATEGY],
        output_path = haskey(options, :output_path) ? abspath(options[:output_path]) : nothing,
    )

    return validate_config(config)
end

strategy_specs() = Common.strategy_specs(String[DEFAULT_MODULI_STRATEGY])
strategy_specs(config::CliConfig) = Common.strategy_specs(config.moduli_strategies)
expression_terms(qe::PC.QuadExpr) = Common.expression_terms(qe)
exact_bound_tightening(con::PC.Constraint, var_bounds::Dict{PC.VarId, PC.IntVar}) =
    Common.exact_bound_tightening(con, var_bounds)

function _sample_nonzero_coefficient(rng::AbstractRNG, coefficients::AbstractVector{Int})
    return rand(rng, coefficients)
end

function random_expression_terms(rng::AbstractRNG, config::CliConfig, nvars::Int)
    coefficients = coefficient_palette(rng, config)
    quad_terms = QuadTerm[]
    lin_terms = LinTerm[]

    sizehint!(quad_terms, nvars * (nvars + 1) ÷ 2)
    sizehint!(lin_terms, nvars)

    for first_id in 1:nvars
        for second_id in first_id:nvars
            rand(rng) < config.density || continue
            coefficient = _sample_nonzero_coefficient(rng, coefficients)
            push!(quad_terms, (Float64(coefficient), first_id, second_id))
        end
    end

    for var_id in 1:nvars
        rand(rng) < config.density || continue
        coefficient = _sample_nonzero_coefficient(rng, coefficients)
        push!(lin_terms, (Float64(coefficient), var_id))
    end

    return quad_terms, lin_terms
end

function _generate_constraint_sample_once(
        rng::AbstractRNG,
        config::CliConfig,
        nvars::Int,
        seed::Int,
        con_id::Int,
        generation_try::Int,
    )
    quad_terms, lin_terms = random_expression_terms(rng, config, nvars)
    isempty(quad_terms) && isempty(lin_terms) && return nothing

    qe = PC.QuadExpr(quad_terms, lin_terms)
    x_star = [rand(rng, config.domain_lb:config.domain_ub) for _ in 1:nvars]
    rhs = PC.eval_full(qe, x_star)
    delta_1 = rand(rng, config.offset_lb:config.offset_ub)
    delta_2 = rand(rng, config.offset_lb:config.offset_ub)
    con = PC.Constraint(con_id, qe, rhs - delta_1, rhs + delta_2)
    model = Common.build_one_constraint_model(nvars, con, config.domain_lb, config.domain_ub)

    PC.normalize!(model; scale_gcd = true)
    (model.infeasible || length(model.cons) != 1) && return nothing

    return ConstraintSample(seed, nvars, model, only(model.cons), x_star, generation_try)
end

function generate_constraint_sample(
        rng::AbstractRNG,
        config::CliConfig,
        nvars::Int;
        seed::Int = 0,
        con_id::Int = 1,
    )
    for generation_try in 1:config.max_generation_tries
        sample = _generate_constraint_sample_once(
            rng,
            config,
            nvars,
            seed,
            con_id,
            generation_try,
        )
        sample === nothing || return sample
    end

    error(
        "failed to generate a nondegenerate quadratic constraint after " *
        "$(config.max_generation_tries) attempts; increase density or max_generation_tries",
    )
end

function generate_constraint_sample(config::CliConfig, nvars::Int, seed::Int; con_id::Int = 1)
    return generate_constraint_sample(
        MersenneTwister(seed),
        config,
        nvars;
        seed = seed,
        con_id = con_id,
    )
end

function constraint_seed(config::CliConfig, n_index::Int, constraint_index::Int)
    return config.seed_base + ((n_index - 1) * config.count + constraint_index) * config.seed_step
end

function result_row(result::StrategyResult)
    return (
        nvars = result.nvars,
        density = result.density,
        domain_lb = result.domain_lb,
        domain_ub = result.domain_ub,
        max_distinct_coeffs = result.max_distinct_coeffs,
        exact_enumeration = result.exact_enumeration,
        Common.result_metric_tail(result)...,
    )
end

function result_rows(results::Vector{StrategyResult})
    rows = NamedTuple[]
    for result in results
        push!(rows, result_row(result))
    end
    return rows
end

function run_experiment(config::CliConfig)
    strategies = strategy_specs(config)
    results = StrategyResult[]
    generated_constraints = Dict{Int, Int}(nvars => 0 for nvars in config.nvars)

    for (n_index, nvars) in enumerate(config.nvars)
        n_results = [
            StrategyResult(
                nvars = nvars,
                density = config.density,
                domain_lb = config.domain_lb,
                domain_ub = config.domain_ub,
                max_distinct_coeffs = config.max_distinct_coeffs,
                name = strategy.name,
                moduli = copy(strategy.moduli),
                exact_enumeration = config.exact_enumeration,
            )
            for strategy in strategies
        ]

        for constraint_index in 0:(config.count - 1)
            seed = constraint_seed(config, n_index, constraint_index)
            sample = generate_constraint_sample(
                config,
                nvars,
                seed;
                con_id = constraint_index + 1,
            )
            generated_constraints[nvars] += 1

            exact = config.exact_enumeration ?
                exact_bound_tightening(sample.con, sample.model.vars) :
                nothing

            for result in n_results
                Common.record_strategy_trial!(
                    result,
                    sample.model,
                    sample.con,
                    exact,
                    config.treewidth_threshold,
                )
            end
        end

        append!(results, n_results)
    end

    rows = result_rows(results)
    return (
        config = config,
        strategies = strategies,
        results = results,
        rows = rows,
        generated_constraints = generated_constraints,
    )
end

write_csv(path::AbstractString, rows) = Common.write_csv(path, rows)

function print_config(result)
    config = result.config
    println("Residue random quadratic constraint experiment")
    println("count = $(config.count)")
    println("nvars = $(join(config.nvars, ","))")
    println("seed_base = $(config.seed_base)")
    println("seed_step = $(config.seed_step)")
    println("density = $(config.density)")
    println("domain = $(config.domain_lb):$(config.domain_ub)")
    println("coeff_range = $(config.coeff_lb):$(config.coeff_ub) excluding 0")
    println("max_distinct_coeffs = $(config.max_distinct_coeffs)")
    println("offset_range = $(config.offset_lb):$(config.offset_ub)")
    println("treewidth_threshold = $(config.treewidth_threshold)")
    println("max_generation_tries = $(config.max_generation_tries)")
    println("moduli_strategies = $(join(config.moduli_strategies, ","))")
    println("exact_enumeration = $(config.exact_enumeration)")
    for nvars in config.nvars
        println("generated_constraints[$nvars] = $(result.generated_constraints[nvars])")
    end
    println()
    return nothing
end

function print_table(rows)
    @printf(
        "%8s %-16s %10s %12s %14s %14s %14s %14s %14s %14s\n",
        "nvars",
        "strategy",
        "n_moduli",
        "constraints",
        "pct_eq",
        "avg_range_red",
        "pct_optimal",
        "pct_tight",
        "avg_gap",
        "avg_time_sec",
    )
    for row in rows
        @printf(
            "%8d %-16s %10d %12d %14.6f %14.6f %14s %14.6f %14s %14.6f\n",
            row.nvars,
            row.strategy,
            row.num_moduli,
            row.constraints,
            row.pct_constraints_tightened_to_equality,
            row.avg_relative_bound_range_reduction,
            Common.metric_float(row.pct_bounds_fully_tightened_to_optimal),
            row.pct_bounds_tightened,
            Common.metric_float(row.avg_bound_gap_to_optimal),
            row.avg_wall_time_sec_per_constraint,
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
        write_csv(config.output_path, result.rows)
        println()
        println("Wrote CSV to $(config.output_path)")
    end

    return result
end

end # module

if _RUNNING_AS_SCRIPT
    ResidueRandomQuadraticConstraintExperiment.main(copy(ARGS))
end
