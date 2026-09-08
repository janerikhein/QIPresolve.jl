#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module ResidueRandomDistanceConstraintExperiment

using Printf
using Random

include("residue_constraint_experiment_common.jl")

const Common = ResidueConstraintExperimentCommon
const PC = Common.PC
const QuadTerm = Common.QuadTerm
const LinTerm = Common.LinTerm

const DEFAULT_COUNT = 1000
const DEFAULT_R = 50
const DEFAULT_ALPHA = 0.01
const DEFAULT_SEED_BASE = 32000
const DEFAULT_SEED_STEP = 1
const DEFAULT_MAX_GENERATION_TRIES = 100
const DEFAULT_EXACT_ENUMERATION = false
const DEFAULT_MODULI_STRATEGIES = Common.moduli_family_strategy_names()
const NVARS = 4

const CLI_KEYS = Dict(
    "count" => :count,
    "r" => :R,
    "alpha" => :alpha,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
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
    R::Int = DEFAULT_R
    alpha::Float64 = DEFAULT_ALPHA
    seed_base::Int = DEFAULT_SEED_BASE
    seed_step::Int = DEFAULT_SEED_STEP
    treewidth_threshold::Int = Common.DEFAULT_TREEWIDTH_THRESHOLD
    max_generation_tries::Int = DEFAULT_MAX_GENERATION_TRIES
    exact_enumeration::Bool = DEFAULT_EXACT_ENUMERATION
    moduli_strategies::Vector{String} = copy(DEFAULT_MODULI_STRATEGIES)
    output_path::Union{Nothing, String} = nothing
end

struct GridPoint
    x::Int
    y::Int
end

struct ConstraintSample
    seed::Int
    R::Int
    alpha::Float64
    model::PC.QPModel
    con::PC.Constraint
    point_1::GridPoint
    point_2::GridPoint
    squared_distance::Int
    x_star::Vector{Int}
    generation_tries::Int
end

Base.@kwdef mutable struct StrategyResult
    R::Int
    alpha::Float64
    nvars::Int
    domain_lb::Int
    domain_ub::Int
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
      julia --project=. scripts/residue_random_distance_constraint_experiment.jl [options]

    Options:
      --count n                       Constraints to generate, default $DEFAULT_COUNT
      --R n                           Coordinate box [-R, R]^2, default $DEFAULT_R
      --alpha p                       Relative squared-distance tolerance, default $DEFAULT_ALPHA
      --seed-base n                   First random seed, default $DEFAULT_SEED_BASE
      --seed-step n                   Seed increment, default $DEFAULT_SEED_STEP
      --treewidth-threshold n         Residue DP treewidth threshold
      --max-generation-tries n        Attempts to avoid degenerate point pairs, default $DEFAULT_MAX_GENERATION_TRIES
      --exact-enumeration bool        Run brute-force exact bound comparison, default $DEFAULT_EXACT_ENUMERATION
      --moduli-strategy list          Comma-separated strategies, default $(join(DEFAULT_MODULI_STRATEGIES, ","))
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

function validate_config(config::CliConfig)
    config.count >= 1 || error("count must be >= 1")
    config.R >= 1 || error("R must be >= 1")
    Common.validate_nonnegative_alpha(config.alpha)
    config.seed_step >= 0 || error("seed_step must be >= 0")
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
        R = Common.parse_int(get(options, :R, string(DEFAULT_R)), "R"),
        alpha = Common.parse_float(get(options, :alpha, string(DEFAULT_ALPHA)), "alpha"),
        seed_base = Common.parse_int(get(options, :seed_base, string(DEFAULT_SEED_BASE)), "seed_base"),
        seed_step = Common.parse_int(get(options, :seed_step, string(DEFAULT_SEED_STEP)), "seed_step"),
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
            copy(DEFAULT_MODULI_STRATEGIES),
        output_path = haskey(options, :output_path) ? abspath(options[:output_path]) : nothing,
    )

    return validate_config(config)
end

strategy_specs() = Common.strategy_specs(DEFAULT_MODULI_STRATEGIES)
strategy_specs(config::CliConfig) = Common.strategy_specs(config.moduli_strategies)
expression_terms(qe::PC.QuadExpr) = Common.expression_terms(qe)
exact_bound_tightening(con::PC.Constraint, var_bounds::Dict{PC.VarId, PC.IntVar}) =
    Common.exact_bound_tightening(con, var_bounds)

@inline function rand_point(rng::AbstractRNG, R::Int)
    return GridPoint(rand(rng, -R:R), rand(rng, -R:R))
end

@inline squared_dist2(p1::GridPoint, p2::GridPoint) = (p1.x - p2.x)^2 + (p1.y - p2.y)^2

function random_distinct_point_pair(rng::AbstractRNG, R::Int, max_generation_tries::Int)
    for generation_try in 1:max_generation_tries
        point_1 = rand_point(rng, R)
        point_2 = rand_point(rng, R)
        point_1 == point_2 && continue
        return point_1, point_2, generation_try
    end

    error(
        "failed to sample two distinct points after $max_generation_tries attempts; " *
        "increase R or max_generation_tries",
    )
end

function distance_expression()
    quad_terms = QuadTerm[
        (1.0, 1, 1),
        (1.0, 2, 2),
        (1.0, 3, 3),
        (1.0, 4, 4),
        (-2.0, 1, 3),
        (-2.0, 2, 4),
    ]
    return PC.QuadExpr(quad_terms, LinTerm[])
end

function generate_constraint_sample(
        rng::AbstractRNG,
        config::CliConfig;
        seed::Int = 0,
        con_id::Int = 1,
    )
    point_1, point_2, generation_tries = random_distinct_point_pair(
        rng,
        config.R,
        config.max_generation_tries,
    )
    squared_distance = squared_dist2(point_1, point_2)
    @assert squared_distance > 0

    lhs = (1.0 - config.alpha) * squared_distance
    rhs = (1.0 + config.alpha) * squared_distance
    con = PC.Constraint(con_id, distance_expression(), lhs, rhs)
    model = Common.build_one_constraint_model(NVARS, con, -config.R, config.R)

    PC.normalize!(model; scale_gcd = true)
    model.infeasible && error("generated distance constraint became infeasible during normalization")
    length(model.cons) == 1 || error("generated distance constraint was eliminated during normalization")

    x_star = [point_1.x, point_1.y, point_2.x, point_2.y]
    return ConstraintSample(
        seed,
        config.R,
        config.alpha,
        model,
        only(model.cons),
        point_1,
        point_2,
        squared_distance,
        x_star,
        generation_tries,
    )
end

function generate_constraint_sample(config::CliConfig, seed::Int; con_id::Int = 1)
    return generate_constraint_sample(
        MersenneTwister(seed),
        config;
        seed = seed,
        con_id = con_id,
    )
end

function constraint_seed(config::CliConfig, constraint_index::Int)
    return config.seed_base + constraint_index * config.seed_step
end

function result_row(result::StrategyResult)
    return (
        R = result.R,
        alpha = result.alpha,
        nvars = result.nvars,
        domain_lb = result.domain_lb,
        domain_ub = result.domain_ub,
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
    results = [
        StrategyResult(
            R = config.R,
            alpha = config.alpha,
            nvars = NVARS,
            domain_lb = -config.R,
            domain_ub = config.R,
            name = strategy.name,
            moduli = copy(strategy.moduli),
            exact_enumeration = config.exact_enumeration,
        )
        for strategy in strategies
    ]
    generated_constraints = 0

    for constraint_index in 0:(config.count - 1)
        seed = constraint_seed(config, constraint_index)
        sample = generate_constraint_sample(
            config,
            seed;
            con_id = constraint_index + 1,
        )
        generated_constraints += 1

        exact = config.exact_enumeration ?
            exact_bound_tightening(sample.con, sample.model.vars) :
            nothing

        for result in results
            Common.record_strategy_trial!(
                result,
                sample.model,
                sample.con,
                exact,
                config.treewidth_threshold,
            )
        end
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
    println("Residue random distance constraint experiment")
    println("count = $(config.count)")
    println("R = $(config.R)")
    println("alpha = $(config.alpha)")
    println("domain = $(-config.R):$(config.R)")
    println("seed_base = $(config.seed_base)")
    println("seed_step = $(config.seed_step)")
    println("treewidth_threshold = $(config.treewidth_threshold)")
    println("max_generation_tries = $(config.max_generation_tries)")
    println("moduli_strategies = $(join(config.moduli_strategies, ","))")
    println("exact_enumeration = $(config.exact_enumeration)")
    println("generated_constraints = $(result.generated_constraints)")
    println()
    return nothing
end

function print_table(rows)
    @printf(
        "%8s %10s %-16s %10s %12s %14s %14s %14s %14s %14s %14s\n",
        "R",
        "alpha",
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
            "%8d %10.6f %-16s %10d %12d %14.6f %14.6f %14s %14.6f %14s %14.6f\n",
            row.R,
            row.alpha,
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
    ResidueRandomDistanceConstraintExperiment.main(copy(ARGS))
end
