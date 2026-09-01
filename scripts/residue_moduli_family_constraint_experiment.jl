#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module ResidueModuliFamilyConstraintExperiment

using CSV
using Printf
using Random

import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC

const DEFAULT_COUNT = 1000
const DEFAULT_NVARS = (3,4,5,6,7,8,9,10)
const DEFAULT_SEED_BASE = 30000
const DEFAULT_SEED_STEP = 1
const DEFAULT_DOMAIN_LB = 0
const DEFAULT_DOMAIN_UB = 1
const DEFAULT_EXTRA_EDGE_PROBABILITY = 0.1
const DEFAULT_COEFF_LB = -50
const DEFAULT_COEFF_UB = 50
const DEFAULT_MAX_DISTINCT_COEFFS = 10
const DEFAULT_DIAG_PROBABILITY = 0.1
const DEFAULT_LINEAR_PROBABILITY = 0.1
const DEFAULT_OFFSET_LB = 1
const DEFAULT_OFFSET_UB = 10
const DEFAULT_MODULUS_LIMIT = 64

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}

const CLI_KEYS = Dict(
    "count" => :count,
    "nvars" => :nvars,
    "n-vars" => :nvars,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "diag-probability" => :diag_probability,
    "diagonal-probability" => :diag_probability,
    "linear-probability" => :linear_probability,
    "lin-probability" => :linear_probability,
    "domain-lb" => :domain_lb,
    "domain-ub" => :domain_ub,
    "extra-edge-probability" => :extra_edge_probability,
    "extra-edge-prob" => :extra_edge_probability,
    "extra-edges" => :extra_edge_probability,
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
    "modulus-limit" => :modulus_limit,
    "moduli-limit" => :modulus_limit,
    "max-modulus" => :modulus_limit,
    "treewidth-threshold" => :treewidth_threshold,
    "exact" => :exact,
    "output" => :output_path,
)

const BOOL_KEYS = Set([:exact])

Base.@kwdef struct CliConfig
    count::Int = DEFAULT_COUNT
    nvars::Vector{Int} = collect(DEFAULT_NVARS)
    seed_base::Int = DEFAULT_SEED_BASE
    seed_step::Int = DEFAULT_SEED_STEP
    diag_probability::Float64 = DEFAULT_DIAG_PROBABILITY
    linear_probability::Float64 = DEFAULT_LINEAR_PROBABILITY
    domain_lb::Int = DEFAULT_DOMAIN_LB
    domain_ub::Int = DEFAULT_DOMAIN_UB
    extra_edge_probability::Float64 = DEFAULT_EXTRA_EDGE_PROBABILITY
    coeff_lb::Int = DEFAULT_COEFF_LB
    coeff_ub::Int = DEFAULT_COEFF_UB
    max_distinct_coeffs::Union{Nothing, Int} = DEFAULT_MAX_DISTINCT_COEFFS
    offset_lb::Int = DEFAULT_OFFSET_LB
    offset_ub::Int = DEFAULT_OFFSET_UB
    modulus_limit::Int = DEFAULT_MODULUS_LIMIT
    treewidth_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD
    exact::Bool = false
    output_path::Union{Nothing, String} = nothing
end

struct StrategySpec
    name::String
    moduli::Vector{Int}
end

struct ConstraintSample
    seed::Int
    nvars::Int
    model::PC.QPModel
    con::PC.Constraint
    x_star::Vector{Int}
end

struct ExactBoundTightening
    lhs::Float64
    rhs::Float64
    relative_bound_reduction::Float64
    assignment_count::Int
end

struct ExpressionTerms
    lin_terms::Vector{Tuple{PC.VarId, Int}}
    quad_terms::Vector{Tuple{PC.VarId, PC.VarId, Int}}
end

Base.@kwdef mutable struct StrategyResult
    nvars::Int
    modulus_limit::Int
    name::String
    moduli::Vector{Int}
    constraints::Int = 0
    bounds_tightened::Int = 0
    total_relative_bound_reduction::Float64 = 0.0
    total_residue_time_sec::Float64 = 0.0
    exact_assignments_per_constraint::Union{Missing, Int} = missing
    total_optimal_relative_bound_reduction::Float64 = 0.0
    bounds_tightened_to_optimal::Int = 0
end

function usage()
    return """
    Usage:
      julia --project=. scripts/residue_moduli_family_constraint_experiment.jl [options]

    Options:
      --count n                       Constraints per n and strategy, default $DEFAULT_COUNT
      --nvars list                    Comma-separated n values, default $(join(DEFAULT_NVARS, ","))
      --seed-base n                   First random seed, default $DEFAULT_SEED_BASE
      --seed-step n                   Seed increment, default $DEFAULT_SEED_STEP
      --diag-probability p            Diagonal x_i^2 term probability, default $DEFAULT_DIAG_PROBABILITY
      --linear-probability p          Linear x_i term probability, default $DEFAULT_LINEAR_PROBABILITY
      --domain-lb n                   Variable lower bound, default $DEFAULT_DOMAIN_LB
      --domain-ub n                   Variable upper bound, default $DEFAULT_DOMAIN_UB
      --extra-edge-probability p      Probability for each non-tree edge, default $DEFAULT_EXTRA_EDGE_PROBABILITY
      --coeff-lb n                    Coefficient lower bound, default $DEFAULT_COEFF_LB
      --coeff-ub n                    Coefficient upper bound, default $DEFAULT_COEFF_UB
      --max-distinct-coeffs n         Max distinct coefficients per generated constraint, default all
      --offset-lb n                   Lower offset bound for sampled constraint slack, default $DEFAULT_OFFSET_LB
      --offset-ub n                   Upper offset bound for sampled constraint slack, default $DEFAULT_OFFSET_UB
      --modulus-limit n               Exclusive upper limit for generated moduli, default $DEFAULT_MODULUS_LIMIT
      --treewidth-threshold n         Residue DP treewidth threshold
      --exact[=true|false]            Compute exact optimal tightening for all nvars values, default false
      --output path                   Optional CSV output path
      -h, --help                      Show this help
    """
end

function parse_int(value::AbstractString, name::AbstractString)::Int
    parsed = tryparse(Int, strip(value))
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function parse_float(value::AbstractString, name::AbstractString)::Float64
    parsed = tryparse(Float64, strip(value))
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function parse_bool(value::AbstractString, name::AbstractString)::Bool
    normalized = lowercase(strip(value))
    normalized in ("true", "1", "yes", "y", "on") && return true
    normalized in ("false", "0", "no", "n", "off") && return false
    error("Invalid $name: $value. Expected true/false, 1/0, yes/no, or on/off.")
end

function parse_int_list(value::AbstractString, name::AbstractString)::Vector{Int}
    normalized = replace(strip(value), " " => "")
    values = Int[]

    for part in split(normalized, ","; keepempty = false)
        if occursin(":", part)
            endpoints = split(part, ":"; limit = 2)
            length(endpoints) == 2 || error("Invalid $name range: $part")
            lower = parse_int(endpoints[1], name)
            upper = parse_int(endpoints[2], name)
            lower <= upper || error("Invalid $name range: lower bound must be <= upper bound")
            append!(values, lower:upper)
        else
            push!(values, parse_int(part, name))
        end
    end

    isempty(values) && error("$name must contain at least one integer")
    return values
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
            key = _lookup_option_key(raw_key)
            if key in BOOL_KEYS && (index == length(args) || startswith(args[index + 1], "--"))
                value = "true"
                options[key] = value
                index += consumed
                continue
            end

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

function validate_probability(name::AbstractString, value::Float64)
    isfinite(value) && 0.0 <= value <= 1.0 ||
        error("$name must be in [0, 1], got $value")
    return value
end

function coefficient_values(config::CliConfig)
    values = [value for value in config.coeff_lb:config.coeff_ub if value != 0]
    isempty(values) && error("coefficient range must contain a nonzero value")
    return values
end

function validate_config(config::CliConfig)
    config.count >= 1 || error("count must be >= 1")
    isempty(config.nvars) && error("nvars must contain at least one value")
    all(>=(1), config.nvars) || error("all nvars values must be >= 1")
    config.seed_step >= 0 || error("seed_step must be >= 0")
    validate_probability("diag_probability", config.diag_probability)
    validate_probability("linear_probability", config.linear_probability)
    config.domain_lb <= config.domain_ub || error("domain_lb must be <= domain_ub")
    validate_probability("extra_edge_probability", config.extra_edge_probability)
    config.coeff_lb <= config.coeff_ub || error("coeff_lb must be <= coeff_ub")
    coefficient_values(config)
    if config.max_distinct_coeffs !== nothing
        config.max_distinct_coeffs >= 1 ||
            error("max_distinct_coeffs must be >= 1")
    end
    config.offset_lb >= 0 || error("offset_lb must be >= 0")
    config.offset_lb <= config.offset_ub || error("offset_lb must be <= offset_ub")
    config.modulus_limit >= 3 || error("modulus_limit must be >= 3")
    config.treewidth_threshold >= 0 || error("treewidth_threshold must be >= 0")

    return config
end

function build_config(args::Vector{String})::Union{Nothing, CliConfig}
    options = parse_raw_options(args)
    options === nothing && return nothing

    config = CliConfig(
        count = parse_int(get(options, :count, string(DEFAULT_COUNT)), "count"),
        nvars = haskey(options, :nvars) ?
            parse_int_list(options[:nvars], "nvars") :
            collect(DEFAULT_NVARS),
        seed_base = parse_int(get(options, :seed_base, string(DEFAULT_SEED_BASE)), "seed_base"),
        seed_step = parse_int(get(options, :seed_step, string(DEFAULT_SEED_STEP)), "seed_step"),
        diag_probability = parse_float(
            get(options, :diag_probability, string(DEFAULT_DIAG_PROBABILITY)),
            "diag_probability",
        ),
        linear_probability = parse_float(
            get(options, :linear_probability, string(DEFAULT_LINEAR_PROBABILITY)),
            "linear_probability",
        ),
        domain_lb = parse_int(get(options, :domain_lb, string(DEFAULT_DOMAIN_LB)), "domain_lb"),
        domain_ub = parse_int(get(options, :domain_ub, string(DEFAULT_DOMAIN_UB)), "domain_ub"),
        extra_edge_probability = parse_float(
            get(options, :extra_edge_probability, string(DEFAULT_EXTRA_EDGE_PROBABILITY)),
            "extra_edge_probability",
        ),
        coeff_lb = parse_int(get(options, :coeff_lb, string(DEFAULT_COEFF_LB)), "coeff_lb"),
        coeff_ub = parse_int(get(options, :coeff_ub, string(DEFAULT_COEFF_UB)), "coeff_ub"),
        max_distinct_coeffs = haskey(options, :max_distinct_coeffs) ?
            parse_int(options[:max_distinct_coeffs], "max_distinct_coeffs") :
            DEFAULT_MAX_DISTINCT_COEFFS,
        offset_lb = parse_int(get(options, :offset_lb, string(DEFAULT_OFFSET_LB)), "offset_lb"),
        offset_ub = parse_int(get(options, :offset_ub, string(DEFAULT_OFFSET_UB)), "offset_ub"),
        modulus_limit = parse_int(
            get(options, :modulus_limit, string(DEFAULT_MODULUS_LIMIT)),
            "modulus_limit",
        ),
        treewidth_threshold = parse_int(
            get(
                options,
                :treewidth_threshold,
                string(QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD),
            ),
            "treewidth_threshold",
        ),
        exact = parse_bool(get(options, :exact, "false"), "exact"),
        output_path = haskey(options, :output_path) ? abspath(options[:output_path]) : nothing,
    )

    return validate_config(config)
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

function primes_less_than(limit::Int)
    return [candidate for candidate in 2:(limit - 1) if is_prime(candidate)]
end

function prime_power_moduli_less_than(limit::Int)
    values = Set{Int}()
    for prime in primes_less_than(limit)
        value = prime
        while value < limit
            push!(values, value)
            value > div(limit - 1, prime) && break
            value *= prime
        end
    end
    return sort!(collect(values))
end

function is_divisibility_antichain(values::AbstractVector{<:Integer})
    for (index, left) in enumerate(values)
        for right in @view values[(index + 1):end]
            (left % right == 0 || right % left == 0) && return false
        end
    end
    return true
end

function full_moduli_less_than(limit::Int)
    upper = limit - 1
    return 2 <= upper ? collect(2:upper) : Int[]
end

function strategy_specs(modulus_limit::Int = DEFAULT_MODULUS_LIMIT)
    return [
        StrategySpec("primes", primes_less_than(modulus_limit)),
        StrategySpec("prime_powers", prime_power_moduli_less_than(modulus_limit)),
        StrategySpec("full", full_moduli_less_than(modulus_limit)),
    ]
end

@inline edge_pair(u::Int, v::Int) = u < v ? (u, v) : (v, u)

function random_spanning_tree_edges(rng::AbstractRNG, nvars::Int)
    nvars >= 1 || throw(ArgumentError("nvars must be >= 1, got $nvars"))
    nvars == 1 && return Tuple{Int, Int}[]

    prufer = [rand(rng, 1:nvars) for _ in 1:max(nvars - 2, 0)]
    degree = ones(Int, nvars)
    for var_id in prufer
        degree[var_id] += 1
    end

    edges = Tuple{Int, Int}[]
    sizehint!(edges, nvars - 1)

    for var_id in prufer
        leaf = findfirst(==(1), degree)
        leaf === nothing && error("failed to decode Prufer sequence")
        push!(edges, edge_pair(leaf, var_id))
        degree[leaf] -= 1
        degree[var_id] -= 1
    end

    leaves = findall(==(1), degree)
    length(leaves) == 2 || error("failed to finish Prufer tree decoding")
    push!(edges, edge_pair(leaves[1], leaves[2]))
    return edges
end

function random_tree_plus_edges(rng::AbstractRNG, nvars::Int, extra_edge_probability::Float64)
    validate_probability("extra_edge_probability", extra_edge_probability)
    tree_edges = random_spanning_tree_edges(rng, nvars)
    edge_set = Set(tree_edges)

    candidates = Tuple{Int, Int}[]
    for u in 1:(nvars - 1)
        for v in (u + 1):nvars
            edge = (u, v)
            edge in edge_set && continue
            push!(candidates, edge)
        end
    end

    edges = copy(tree_edges)
    for edge in candidates
        rand(rng) < extra_edge_probability || continue
        push!(edges, edge)
    end

    return edges
end

function _sample_nonzero_coefficient(rng::AbstractRNG, coefficients::AbstractVector{Int})
    return rand(rng, coefficients)
end

function coefficient_palette(rng::AbstractRNG, config::CliConfig)
    coefficients = coefficient_values(config)
    config.max_distinct_coeffs === nothing && return coefficients

    palette_size = min(config.max_distinct_coeffs, length(coefficients))
    palette_size == length(coefficients) && return coefficients

    order = randperm(rng, length(coefficients))
    palette = coefficients[order[1:palette_size]]
    return sort!(palette)
end

function build_one_constraint_model(
        nvars::Int,
        con::PC.Constraint,
        domain_lb::Int,
        domain_ub::Int,
    )
    vars = Dict{PC.VarId, PC.IntVar}()
    for var_id in 1:nvars
        vars[var_id] = PC.IntVar(Float64(domain_lb), Float64(domain_ub))
    end
    return PC.QPModel(vars, [con], PC.QuadExpr(QuadTerm[], LinTerm[]), :min)
end

function generate_constraint_sample(
        rng::AbstractRNG,
        config::CliConfig,
        nvars::Int;
        seed::Int = 0,
        con_id::Int = 1,
    )
    coefficients = coefficient_palette(rng, config)
    edges = random_tree_plus_edges(rng, nvars, config.extra_edge_probability)

    quad_terms = QuadTerm[]
    sizehint!(quad_terms, length(edges) + nvars)
    for (u, v) in edges
        coefficient = _sample_nonzero_coefficient(rng, coefficients)
        push!(quad_terms, (Float64(coefficient), u, v))
    end

    for var_id in 1:nvars
        rand(rng) < config.diag_probability || continue
        coefficient = _sample_nonzero_coefficient(rng, coefficients)
        push!(quad_terms, (Float64(coefficient), var_id, var_id))
    end

    lin_terms = LinTerm[]
    sizehint!(lin_terms, nvars)
    for var_id in 1:nvars
        rand(rng) < config.linear_probability || continue
        coefficient = _sample_nonzero_coefficient(rng, coefficients)
        push!(lin_terms, (Float64(coefficient), var_id))
    end

    qe = PC.QuadExpr(quad_terms, lin_terms)
    x_star = [rand(rng, config.domain_lb:config.domain_ub) for _ in 1:nvars]
    rhs = PC.eval_full(qe, x_star)
    delta_1 = rand(rng, config.offset_lb:config.offset_ub)
    delta_2 = rand(rng, config.offset_lb:config.offset_ub)
    con = PC.Constraint(con_id, qe, rhs - delta_1, rhs + delta_2)
    model = build_one_constraint_model(nvars, con, config.domain_lb, config.domain_ub)

    PC.normalize!(model; scale_gcd = true)
    model.infeasible && error("generated constraint became infeasible during normalization")
    length(model.cons) == 1 || error("generated constraint was eliminated during normalization")

    return ConstraintSample(seed, nvars, model, only(model.cons), x_star)
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

function one_constraint_model(model::PC.QPModel, con::PC.Constraint)
    return PC.QPModel(model.vars, [con], model.obj_expr, model.obj_sense)
end

bound_snapshot(con::PC.Constraint) = (lhs = con.lhs, rhs = con.rhs)

function relative_bound_reduction(before, lhs::Real, rhs::Real)
    isfinite(before.lhs) && isfinite(before.rhs) || return 0.0

    baseline_range = before.rhs - before.lhs
    baseline_range > 0.0 || return 0.0

    improvement = 0.0
    isfinite(lhs) && (improvement += max(0.0, Float64(lhs) - before.lhs))
    isfinite(rhs) && (improvement += max(0.0, before.rhs - Float64(rhs)))
    return improvement / baseline_range
end

relative_bound_reduction(before, con::PC.Constraint) =
    relative_bound_reduction(before, con.lhs, con.rhs)

function tightened_bound_count(before, con::PC.Constraint)
    count = 0
    isfinite(before.lhs) && con.lhs > before.lhs && (count += 1)
    isfinite(before.rhs) && con.rhs < before.rhs && (count += 1)
    return count
end

function run_residue_strategy!(
        model::PC.QPModel,
        con::PC.Constraint,
        moduli::AbstractVector{<:Integer},
        treewidth_threshold::Int,
    )
    model.infeasible && return true

    status = PC._residue_constraint_status(con, model.vars)
    if status == :infeasible
        model.infeasible = true
        return true
    elseif status == :skip || isempty(moduli)
        return false
    end

    standardized = PC._standardize_residue_constraint(con, model.vars, model.obj_expr)
    cache = PC._ResidueCacheEntry[]

    for modulus_value in moduli
        modulus = Int(modulus_value)
        residue_result = PC._compute_achievable_residues(
            modulus,
            standardized.con,
            standardized.var_bounds;
            treewidth_threshold = treewidth_threshold,
        )
        residue_result.saturated && continue

        push!(cache, PC._ResidueCacheEntry(modulus, residue_result.residues))
        PC._reapply_residue_cache!(con, cache, standardized.constraint_shift)

        if con.lhs > con.rhs
            model.infeasible = true
            return true
        end
    end

    return false
end

function _integer_coefficient(value::Float64, description::String)
    isinteger(value) || throw(ArgumentError("$description must be integer-valued, got $value"))
    return trunc(Int, value)
end

function expression_terms(qe::PC.QuadExpr)
    var_ids = sort!(collect(PC.vars(qe)))
    lin_terms = Tuple{PC.VarId, Int}[]
    quad_terms = Tuple{PC.VarId, PC.VarId, Int}[]

    for (index, var_id) in enumerate(var_ids)
        lin_coeff = _integer_coefficient(
            PC.get_lin_coeff(qe, var_id),
            "linear coefficient for $var_id",
        )
        lin_coeff == 0 || push!(lin_terms, (var_id, lin_coeff))

        for other_id in @view var_ids[index:end]
            quad_coeff = _integer_coefficient(
                PC.get_quad_coeff(qe, var_id, other_id),
                "quadratic coefficient for ($var_id, $other_id)",
            )
            quad_coeff == 0 || push!(quad_terms, (var_id, other_id, quad_coeff))
        end
    end

    return ExpressionTerms(lin_terms, quad_terms)
end

function eval_terms(terms::ExpressionTerms, values::AbstractVector{Int})
    total = 0
    @inbounds for (var_id, coefficient) in terms.lin_terms
        total += coefficient * values[var_id]
    end
    @inbounds for (first_id, second_id, coefficient) in terms.quad_terms
        total += coefficient * values[first_id] * values[second_id]
    end
    return total
end

function exact_domains(var_bounds::Dict{PC.VarId, PC.IntVar})
    var_ids = sort!(collect(keys(var_bounds)))
    lbs = Int[]
    ubs = Int[]
    sizehint!(lbs, length(var_ids))
    sizehint!(ubs, length(var_ids))

    for var_id in var_ids
        var = var_bounds[var_id]
        isfinite(var.lb) || throw(ArgumentError("variable $var_id must have finite lower bound"))
        isfinite(var.ub) || throw(ArgumentError("variable $var_id must have finite upper bound"))
        isinteger(var.lb) || throw(ArgumentError("variable $var_id lower bound must be integer-valued"))
        isinteger(var.ub) || throw(ArgumentError("variable $var_id upper bound must be integer-valued"))
        var.lb <= var.ub || throw(ArgumentError("variable $var_id has inconsistent bounds"))
        push!(lbs, trunc(Int, var.lb))
        push!(ubs, trunc(Int, var.ub))
    end

    return var_ids, lbs, ubs
end

function assignment_count(lbs::AbstractVector{Int}, ubs::AbstractVector{Int})
    count = 1
    for (lb, ub) in zip(lbs, ubs)
        count *= ub - lb + 1
    end
    return count
end

function advance_assignment!(
        values::Vector{Int},
        var_ids::AbstractVector{PC.VarId},
        lbs::AbstractVector{Int},
        ubs::AbstractVector{Int},
    )
    for position in length(var_ids):-1:1
        var_id = var_ids[position]
        if values[var_id] < ubs[position]
            values[var_id] += 1
            return true
        end
        values[var_id] = lbs[position]
    end
    return false
end

function exact_bound_tightening(con::PC.Constraint, var_bounds::Dict{PC.VarId, PC.IntVar})
    isfinite(con.lhs) || throw(ArgumentError("constraint lower bound must be finite"))
    isfinite(con.rhs) || throw(ArgumentError("constraint upper bound must be finite"))
    con.lhs <= con.rhs || throw(ArgumentError("constraint has inconsistent bounds"))

    var_ids, lbs, ubs = exact_domains(var_bounds)
    total_assignments = assignment_count(lbs, ubs)
    terms = expression_terms(con.qe)
    values = zeros(Int, maximum(var_ids; init = 0))

    for (position, var_id) in enumerate(var_ids)
        values[var_id] = lbs[position]
    end

    lower_threshold = ceil(Int, con.lhs)
    upper_threshold = floor(Int, con.rhs)
    best_lhs = typemax(Int)
    best_rhs = typemin(Int)

    for _ in 1:total_assignments
        value = eval_terms(terms, values)
        if value >= lower_threshold && value < best_lhs
            best_lhs = value
        end
        if value <= upper_threshold && value > best_rhs
            best_rhs = value
        end
        advance_assignment!(values, var_ids, lbs, ubs)
    end

    best_lhs == typemax(Int) && error("no assignment attains a value >= lower bound")
    best_rhs == typemin(Int) && error("no assignment attains a value <= upper bound")

    before = bound_snapshot(con)
    lhs = Float64(best_lhs)
    rhs = Float64(best_rhs)
    return ExactBoundTightening(
        lhs,
        rhs,
        relative_bound_reduction(before, lhs, rhs),
        total_assignments,
    )
end

function bounds_tightened_to_optimal_count(before, con::PC.Constraint, exact::ExactBoundTightening)
    count = 0
    isfinite(before.lhs) && con.lhs > before.lhs && con.lhs == exact.lhs && (count += 1)
    isfinite(before.rhs) && con.rhs < before.rhs && con.rhs == exact.rhs && (count += 1)
    return count
end

function record_strategy_trial!(
        result::StrategyResult,
        sample::ConstraintSample,
        exact::Union{Nothing, ExactBoundTightening},
        treewidth_threshold::Int,
    )
    trial_con = deepcopy(sample.con)
    trial_model = one_constraint_model(sample.model, trial_con)
    before = bound_snapshot(trial_con)

    start_time = time()
    run_residue_strategy!(
        trial_model,
        trial_con,
        result.moduli,
        treewidth_threshold,
    )
    result.total_residue_time_sec += time() - start_time

    result.constraints += 1
    result.bounds_tightened += tightened_bound_count(before, trial_con)
    result.total_relative_bound_reduction += relative_bound_reduction(before, trial_con)

    if exact !== nothing
        result.exact_assignments_per_constraint = exact.assignment_count
        result.total_optimal_relative_bound_reduction += exact.relative_bound_reduction
        result.bounds_tightened_to_optimal += bounds_tightened_to_optimal_count(before, trial_con, exact)
    end

    return result
end

rate(numerator::Real, denominator::Int) = denominator == 0 ? 0.0 : numerator / denominator
rate(::Missing, ::Int) = missing

function base_row(result::StrategyResult)
    bounds_considered = 2 * result.constraints
    return (
        nvars = result.nvars,
        modulus_limit = result.modulus_limit,
        strategy = result.name,
        moduli = join(result.moduli, " "),
        num_moduli = length(result.moduli),
        constraints = result.constraints,
        bounds_considered = bounds_considered,
        bounds_tightened = result.bounds_tightened,
        fraction_bounds_tightened = rate(result.bounds_tightened, bounds_considered),
        avg_relative_bound_reduction = rate(
            result.total_relative_bound_reduction,
            result.constraints,
        ),
        total_residue_time_sec = result.total_residue_time_sec,
        avg_wall_time_sec_per_constraint = rate(
            result.total_residue_time_sec,
            result.constraints,
        ),
    )
end

function row_with_exact(result::StrategyResult)
    base = base_row(result)
    has_exact = !ismissing(result.exact_assignments_per_constraint)
    bounds_considered = 2 * result.constraints
    return merge(
        base,
        (
            exact_assignments_per_constraint = has_exact ?
                result.exact_assignments_per_constraint :
                missing,
            optimal_avg_relative_bound_reduction = has_exact ?
                rate(result.total_optimal_relative_bound_reduction, result.constraints) :
                missing,
            bounds_tightened_to_optimal = has_exact ?
                result.bounds_tightened_to_optimal :
                missing,
            fraction_bounds_tightened_to_optimal = has_exact ?
                rate(result.bounds_tightened_to_optimal, bounds_considered) :
                missing,
        ),
    )
end

function result_rows(results::Vector{StrategyResult}, include_exact::Bool)
    rows = NamedTuple[]
    for result in results
        push!(rows, include_exact ? row_with_exact(result) : base_row(result))
    end
    return rows
end

function run_experiment(config::CliConfig)
    strategies = strategy_specs(config.modulus_limit)
    results = StrategyResult[]
    generated_constraints = Dict{Int, Int}(nvars => 0 for nvars in config.nvars)

    for (n_index, nvars) in enumerate(config.nvars)
        n_results = [
            StrategyResult(
                nvars = nvars,
                modulus_limit = config.modulus_limit,
                name = strategy.name,
                moduli = copy(strategy.moduli),
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

            exact = config.exact ?
                exact_bound_tightening(sample.con, sample.model.vars) :
                nothing

            for result in n_results
                record_strategy_trial!(
                    result,
                    sample,
                    exact,
                    config.treewidth_threshold,
                )
            end
        end

        append!(results, n_results)
    end

    rows = result_rows(results, config.exact)
    return (
        config = config,
        strategies = strategies,
        results = results,
        rows = rows,
        generated_constraints = generated_constraints,
    )
end

function write_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    CSV.write(path, rows)
    return path
end

function print_config(result)
    config = result.config
    println("Residue moduli family constraint experiment")
    println("count = $(config.count)")
    println("nvars = $(join(config.nvars, ","))")
    println("seed_base = $(config.seed_base)")
    println("seed_step = $(config.seed_step)")
    println("diag_probability = $(config.diag_probability)")
    println("linear_probability = $(config.linear_probability)")
    println("domain = $(config.domain_lb):$(config.domain_ub)")
    println("extra_edge_probability = $(config.extra_edge_probability)")
    println("coeff_range = $(config.coeff_lb):$(config.coeff_ub) excluding 0")
    max_distinct_coeffs = config.max_distinct_coeffs === nothing ?
        "all" :
        string(config.max_distinct_coeffs)
    println("max_distinct_coeffs = $max_distinct_coeffs")
    println("offset_range = $(config.offset_lb):$(config.offset_ub)")
    println("modulus_limit = $(config.modulus_limit)")
    println("treewidth_threshold = $(config.treewidth_threshold)")
    println("exact = $(config.exact)")
    for nvars in config.nvars
        println("generated_constraints[$nvars] = $(result.generated_constraints[nvars])")
    end
    println()
    return nothing
end

function _fmt_missing(value)
    return ismissing(value) ? "" : string(value)
end

function _fmt_float(value)
    return ismissing(value) ? "" : @sprintf("%.6f", value)
end

function print_table(rows; include_exact::Bool)
    if include_exact
        @printf(
            "%8s %-26s %10s %12s %12s %14s %14s %14s %14s %14s\n",
            "nvars",
            "strategy",
            "n_moduli",
            "constraints",
            "bounds_tight",
            "frac_bounds",
            "avg_rel_red",
            "avg_time_sec",
            "opt_avg_red",
            "frac_opt",
        )
        for row in rows
            @printf(
                "%8d %-26s %10d %12d %12d %14.6f %14.6f %14.6f %14s %14s\n",
                row.nvars,
                row.strategy,
                row.num_moduli,
                row.constraints,
                row.bounds_tightened,
                row.fraction_bounds_tightened,
                row.avg_relative_bound_reduction,
                row.avg_wall_time_sec_per_constraint,
                _fmt_float(row.optimal_avg_relative_bound_reduction),
                _fmt_float(row.fraction_bounds_tightened_to_optimal),
            )
        end
    else
        @printf(
            "%8s %-26s %10s %12s %12s %14s %14s %14s\n",
            "nvars",
            "strategy",
            "n_moduli",
            "constraints",
            "bounds_tight",
            "frac_bounds",
            "avg_rel_red",
            "avg_time_sec",
        )
        for row in rows
            @printf(
                "%8d %-26s %10d %12d %12d %14.6f %14.6f %14.6f\n",
                row.nvars,
                row.strategy,
                row.num_moduli,
                row.constraints,
                row.bounds_tightened,
                row.fraction_bounds_tightened,
                row.avg_relative_bound_reduction,
                row.avg_wall_time_sec_per_constraint,
            )
        end
    end
    return nothing
end

function main(args::Vector{String} = copy(ARGS))
    config = build_config(args)
    config === nothing && return nothing

    result = run_experiment(config)
    print_config(result)
    print_table(result.rows; include_exact = config.exact)

    if config.output_path !== nothing
        write_csv(config.output_path, result.rows)
        println()
        println("Wrote CSV to $(config.output_path)")
    end

    return result
end

end # module

if _RUNNING_AS_SCRIPT
    ResidueModuliFamilyConstraintExperiment.main(copy(ARGS))
end
