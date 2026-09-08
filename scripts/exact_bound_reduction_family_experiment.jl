#!/usr/bin/env julia

const _EXACT_BOUND_REDUCTION_FAMILY_RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _EXACT_BOUND_REDUCTION_FAMILY_RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module ExactBoundReductionFamilyExperiment

using CSV
import Graphs
using Printf
using Random
using Statistics

import QIPresolve.PresolvingCore as PC

const DEFAULT_COUNT = 100
const DEFAULT_NVARS = (5,)
const DEFAULT_SEED_BASE = 30000
const DEFAULT_SEED_STEP = 1
const DEFAULT_DOMAIN_LB = 0
const DEFAULT_DOMAIN_UBS = (3,)
const DEFAULT_DENSITY = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
const DEFAULT_COEFF_LB = -50
const DEFAULT_COEFF_UB = 50
const DEFAULT_MAX_DISTINCT_COEFFS = (50,)
const DEFAULT_OFFSET_LB = 1
const DEFAULT_OFFSET_UB = 10

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}

const CLI_KEYS = Dict(
    "count" => :count,
    "nvars" => :nvars,
    "n-vars" => :nvars,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "density" => :densities,
    "densities" => :densities,
    "domain-lb" => :domain_lb,
    "domain-ub" => :domain_ubs,
    "domain-ubs" => :domain_ubs,
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
    "output" => :output_dir,
    "output-dir" => :output_dir,
)

Base.@kwdef struct CliConfig
    count::Int = DEFAULT_COUNT
    nvars::Vector{Int} = collect(DEFAULT_NVARS)
    seed_base::Int = DEFAULT_SEED_BASE
    seed_step::Int = DEFAULT_SEED_STEP
    densities::Vector{Float64} = collect(DEFAULT_DENSITY)
    domain_lb::Int = DEFAULT_DOMAIN_LB
    domain_ubs::Vector{Int} = collect(DEFAULT_DOMAIN_UBS)
    coeff_lb::Int = DEFAULT_COEFF_LB
    coeff_ub::Int = DEFAULT_COEFF_UB
    max_distinct_coeffs::Vector{Int} = collect(DEFAULT_MAX_DISTINCT_COEFFS)
    offset_lb::Int = DEFAULT_OFFSET_LB
    offset_ub::Int = DEFAULT_OFFSET_UB
    output_dir::Union{Nothing, String} = nothing
end

struct ConstraintSample
    seed::Int
    nvars::Int
    domain_ub::Int
    max_distinct_coeffs::Int
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

Base.@kwdef mutable struct SweepResult
    nvars::Int
    density::Float64
    domain_lb::Int
    domain_ub::Int
    max_distinct_coeffs::Int
    constraints::Int = 0
    optimal_relative_bound_reductions::Vector{Float64} = Float64[]
    treewidths::Vector{Int} = Int[]
    total_exact_time_sec::Float64 = 0.0
end

function usage()
    return """
    Usage:
      julia --project=. scripts/exact_bound_reduction_family_experiment.jl [options]

    Options:
      --count n                       Constraints per parameter combination, default $DEFAULT_COUNT
      --nvars list                    Comma-separated n values, default $(join(DEFAULT_NVARS, ","))
      --density, --densities list     Comma-separated density values, default $(join(DEFAULT_DENSITY, ","))
      --domain-ubs list               Comma-separated upper bounds, default $(join(DEFAULT_DOMAIN_UBS, ","))
      --seed-base n                   First random seed, default $DEFAULT_SEED_BASE
      --seed-step n                   Seed increment, default $DEFAULT_SEED_STEP
      --domain-lb n                   Variable lower bound, default $DEFAULT_DOMAIN_LB
      --coeff-lb n                    Coefficient lower bound, default $DEFAULT_COEFF_LB
      --coeff-ub n                    Coefficient upper bound, default $DEFAULT_COEFF_UB
      --max-distinct-coeffs list      Comma-separated max distinct coefficients, default $(join(DEFAULT_MAX_DISTINCT_COEFFS, ","))
      --offset-lb n                   Lower offset bound for sampled constraint slack, default $DEFAULT_OFFSET_LB
      --offset-ub n                   Upper offset bound for sampled constraint slack, default $DEFAULT_OFFSET_UB
      --output, --output-dir path     Optional CSV output directory
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

function parse_float_list(value::AbstractString, name::AbstractString)::Vector{Float64}
    normalized = replace(strip(value), " " => "")
    values = Float64[]

    for part in split(normalized, ","; keepempty = false)
        push!(values, parse_float(part, name))
    end

    isempty(values) && error("$name must contain at least one float")
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
    isempty(config.densities) && error("densities must contain at least one value")
    foreach(density -> validate_probability("density", density), config.densities)
    config.domain_lb <= minimum(config.domain_ubs) ||
        error("all domain_ubs values must be >= domain_lb")
    config.coeff_lb <= config.coeff_ub || error("coeff_lb must be <= coeff_ub")
    coefficient_values(config)
    isempty(config.max_distinct_coeffs) &&
        error("max_distinct_coeffs must contain at least one value")
    all(>=(1), config.max_distinct_coeffs) ||
        error("all max_distinct_coeffs values must be >= 1")
    config.offset_lb >= 0 || error("offset_lb must be >= 0")
    config.offset_lb <= config.offset_ub || error("offset_lb must be <= offset_ub")

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
        densities = haskey(options, :densities) ?
            parse_float_list(options[:densities], "densities") :
            collect(DEFAULT_DENSITY),
        domain_lb = parse_int(get(options, :domain_lb, string(DEFAULT_DOMAIN_LB)), "domain_lb"),
        domain_ubs = haskey(options, :domain_ubs) ?
            parse_int_list(options[:domain_ubs], "domain_ubs") :
            collect(DEFAULT_DOMAIN_UBS),
        coeff_lb = parse_int(get(options, :coeff_lb, string(DEFAULT_COEFF_LB)), "coeff_lb"),
        coeff_ub = parse_int(get(options, :coeff_ub, string(DEFAULT_COEFF_UB)), "coeff_ub"),
        max_distinct_coeffs = haskey(options, :max_distinct_coeffs) ?
            parse_int_list(options[:max_distinct_coeffs], "max_distinct_coeffs") :
            collect(DEFAULT_MAX_DISTINCT_COEFFS),
        offset_lb = parse_int(get(options, :offset_lb, string(DEFAULT_OFFSET_LB)), "offset_lb"),
        offset_ub = parse_int(get(options, :offset_ub, string(DEFAULT_OFFSET_UB)), "offset_ub"),
        output_dir = haskey(options, :output_dir) ? abspath(options[:output_dir]) : nothing,
    )

    return validate_config(config)
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

function random_tree_plus_edges(rng::AbstractRNG, nvars::Int, density::Float64)
    validate_probability("density", density)
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
        rand(rng) < density || continue
        push!(edges, edge)
    end

    return edges
end

function _sample_nonzero_coefficient(rng::AbstractRNG, coefficients::AbstractVector{Int})
    return rand(rng, coefficients)
end

function coefficient_palette(
        rng::AbstractRNG,
        config::CliConfig,
        max_distinct_coeffs::Int,
    )
    coefficients = coefficient_values(config)
    palette_size = min(max_distinct_coeffs, length(coefficients))
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
        nvars::Int,
        domain_ub::Int,
        max_distinct_coeffs::Int;
        density::Float64 = only(config.densities),
        seed::Int = 0,
        con_id::Int = 1,
    )
    validate_probability("density", density)
    coefficients = coefficient_palette(rng, config, max_distinct_coeffs)
    edges = random_tree_plus_edges(rng, nvars, density)

    quad_terms = QuadTerm[]
    sizehint!(quad_terms, length(edges) + nvars)
    for (u, v) in edges
        coefficient = _sample_nonzero_coefficient(rng, coefficients)
        push!(quad_terms, (Float64(coefficient), u, v))
    end

    for var_id in 1:nvars
        rand(rng) < density || continue
        coefficient = _sample_nonzero_coefficient(rng, coefficients)
        push!(quad_terms, (Float64(coefficient), var_id, var_id))
    end

    lin_terms = LinTerm[]
    sizehint!(lin_terms, nvars)
    for var_id in 1:nvars
        rand(rng) < density || continue
        coefficient = _sample_nonzero_coefficient(rng, coefficients)
        push!(lin_terms, (Float64(coefficient), var_id))
    end

    qe = PC.QuadExpr(quad_terms, lin_terms)
    x_star = [rand(rng, config.domain_lb:domain_ub) for _ in 1:nvars]
    rhs = PC.eval_full(qe, x_star)
    delta_1 = rand(rng, config.offset_lb:config.offset_ub)
    delta_2 = rand(rng, config.offset_lb:config.offset_ub)
    con = PC.Constraint(con_id, qe, rhs - delta_1, rhs + delta_2)
    model = build_one_constraint_model(nvars, con, config.domain_lb, domain_ub)

    PC.normalize!(model; scale_gcd = true)
    model.infeasible && error("generated constraint became infeasible during normalization")
    length(model.cons) == 1 || error("generated constraint was eliminated during normalization")

    return ConstraintSample(
        seed,
        nvars,
        domain_ub,
        max_distinct_coeffs,
        model,
        only(model.cons),
        x_star,
    )
end

function generate_constraint_sample(
        config::CliConfig,
        nvars::Int,
        domain_ub::Int,
        max_distinct_coeffs::Int,
        seed::Int;
        density::Float64 = only(config.densities),
        con_id::Int = 1,
    )
    return generate_constraint_sample(
        MersenneTwister(seed),
        config,
        nvars,
        domain_ub,
        max_distinct_coeffs;
        density = density,
        seed = seed,
        con_id = con_id,
    )
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

function combination_seed(config::CliConfig, combination_index::Int, constraint_index::Int)
    return config.seed_base +
        ((combination_index - 1) * config.count + constraint_index) * config.seed_step
end

function _component_graph(graph::Graphs.SimpleGraph{Int}, component::Vector{Int})
    component_pos = Dict(vertex => index for (index, vertex) in enumerate(component))
    component_graph = Graphs.SimpleGraph{Int}(length(component))

    for (first_index, first_vertex) in enumerate(component)
        for second_vertex in @view component[(first_index + 1):end]
            Graphs.has_edge(graph, first_vertex, second_vertex) || continue
            Graphs.add_edge!(
                component_graph,
                component_pos[first_vertex],
                component_pos[second_vertex],
            )
        end
    end

    return component_graph
end

function _treewidth_component(
        graph::Graphs.SimpleGraph{Int},
        component::Vector{Int},
        var_ids::Vector{PC.VarId},
    )
    component_var_ids = var_ids[component]
    var_id_to_pos = Dict(var_id => position for (position, var_id) in enumerate(component_var_ids))
    pc_component = PC.NonSingleton(
        _component_graph(graph, component),
        var_id_to_pos,
        component_var_ids,
        Dict{PC.VarId, Int}(),
        Dict{Tuple{PC.VarId, PC.VarId}, Int}(),
    )
    return PC.tree_decomposition_width(PC.minimum_degree_tree_decomposition(pc_component))
end

function constraint_treewidth(con::PC.Constraint)
    var_ids = sort!(collect(PC.vars(con.qe)))
    length(var_ids) <= 1 && return 0

    graph = Graphs.SimpleGraph{Int}(length(var_ids))
    for first_index in 1:(length(var_ids) - 1)
        first_var_id = var_ids[first_index]
        for second_index in (first_index + 1):length(var_ids)
            second_var_id = var_ids[second_index]
            PC.get_quad_coeff(con.qe, first_var_id, second_var_id) == 0.0 && continue
            Graphs.add_edge!(graph, first_index, second_index)
        end
    end

    Graphs.ne(graph) == 0 && return 0

    width = 0
    for component in Graphs.connected_components(graph)
        length(component) <= 1 && continue
        width = max(width, _treewidth_component(graph, component, var_ids))
    end
    return width
end

function record_trial!(result::SweepResult, sample::ConstraintSample)
    start_time = time()
    exact = exact_bound_tightening(sample.con, sample.model.vars)
    result.total_exact_time_sec += time() - start_time

    result.constraints += 1
    push!(result.optimal_relative_bound_reductions, exact.relative_bound_reduction)
    push!(result.treewidths, constraint_treewidth(sample.con))
    return result
end

rate(numerator::Real, denominator::Int) = denominator == 0 ? 0.0 : numerator / denominator

function treewidth_bucket_percentages(treewidths::AbstractVector{Int}, constraints::Int)
    return (
        tw_1 = rate(100.0 * count(==(1), treewidths), constraints),
        tw_2 = rate(100.0 * count(==(2), treewidths), constraints),
        tw_3 = rate(100.0 * count(==(3), treewidths), constraints),
        tw_ge4 = rate(100.0 * count(>=(4), treewidths), constraints),
    )
end

function result_row(result::SweepResult)
    treewidth_buckets = treewidth_bucket_percentages(result.treewidths, result.constraints)
    return (
        nvars = result.nvars,
        density = result.density,
        domain_lb = result.domain_lb,
        domain_ub = result.domain_ub,
        max_distinct_coeffs = result.max_distinct_coeffs,
        constraints = result.constraints,
        opt_avg_red = mean(result.optimal_relative_bound_reductions),
        opt_std_red = std(result.optimal_relative_bound_reductions),
        treewidth_buckets...,
        avg_wall_time_sec_per_constraint = rate(
            result.total_exact_time_sec,
            result.constraints,
        ),
    )
end

function result_rows(results::Vector{SweepResult})
    rows = NamedTuple[]
    for result in results
        push!(rows, result_row(result))
    end
    return rows
end

function grouped_result_row(
        field::Symbol,
        value,
        constraints::Int,
        reductions,
        treewidths,
        total_time::Float64,
    )
    treewidth_buckets = treewidth_bucket_percentages(treewidths, constraints)
    return NamedTuple{
        (
            field,
            :constraints,
            :opt_avg_red,
            :opt_std_red,
            :tw_1,
            :tw_2,
            :tw_3,
            :tw_ge4,
            :avg_wall_time_sec_per_constraint,
        ),
    }((
        value,
        constraints,
        mean(reductions),
        std(reductions),
        treewidth_buckets.tw_1,
        treewidth_buckets.tw_2,
        treewidth_buckets.tw_3,
        treewidth_buckets.tw_ge4,
        rate(total_time, constraints),
    ))
end

function grouped_result_rows(results::Vector{SweepResult}, field::Symbol)
    reductions_by_value = Dict{Any, Vector{Float64}}()
    treewidths_by_value = Dict{Any, Vector{Int}}()
    constraints_by_value = Dict{Any, Int}()
    total_time_by_value = Dict{Any, Float64}()
    values = Any[]

    for result in results
        value = getproperty(result, field)
        if !haskey(reductions_by_value, value)
            reductions_by_value[value] = Float64[]
            treewidths_by_value[value] = Int[]
            constraints_by_value[value] = 0
            total_time_by_value[value] = 0.0
            push!(values, value)
        end

        append!(reductions_by_value[value], result.optimal_relative_bound_reductions)
        append!(treewidths_by_value[value], result.treewidths)
        constraints_by_value[value] += result.constraints
        total_time_by_value[value] += result.total_exact_time_sec
    end

    rows = NamedTuple[]
    for value in values
        push!(
            rows,
            grouped_result_row(
                field,
                value,
                constraints_by_value[value],
                reductions_by_value[value],
                treewidths_by_value[value],
                total_time_by_value[value],
            ),
        )
    end
    return rows
end

function aggregate_result_rows(results::Vector{SweepResult})
    return (
        by_nvars = grouped_result_rows(results, :nvars),
        by_density = grouped_result_rows(results, :density),
        by_domain_ub = grouped_result_rows(results, :domain_ub),
        by_max_distinct_coeffs = grouped_result_rows(results, :max_distinct_coeffs),
    )
end

function run_experiment(config::CliConfig)
    results = SweepResult[]
    generated_constraints = Dict{Tuple{Int, Float64, Int, Int}, Int}()
    combination_index = 1

    for nvars in config.nvars
        for density in config.densities
            for domain_ub in config.domain_ubs
                for max_distinct_coeffs in config.max_distinct_coeffs
                    result = SweepResult(
                        nvars = nvars,
                        density = density,
                        domain_lb = config.domain_lb,
                        domain_ub = domain_ub,
                        max_distinct_coeffs = max_distinct_coeffs,
                    )
                    key = (nvars, density, domain_ub, max_distinct_coeffs)
                    generated_constraints[key] = 0

                    for constraint_index in 0:(config.count - 1)
                        seed = combination_seed(config, combination_index, constraint_index)
                        sample = generate_constraint_sample(
                            config,
                            nvars,
                            domain_ub,
                            max_distinct_coeffs,
                            seed;
                            density = density,
                            con_id = constraint_index + 1,
                        )
                        generated_constraints[key] += 1
                        record_trial!(result, sample)
                    end

                    push!(results, result)
                    combination_index += 1
                end
            end
        end
    end

    rows = result_rows(results)
    aggregate_rows = aggregate_result_rows(results)
    return (
        config = config,
        results = results,
        rows = rows,
        aggregate_rows = aggregate_rows,
        generated_constraints = generated_constraints,
    )
end

function write_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    CSV.write(path, rows)
    return path
end

function aggregate_output_paths(output_dir::AbstractString)
    return (
        by_nvars = joinpath(output_dir, "exact_res_by_nvars.csv"),
        by_density = joinpath(output_dir, "exact_res_by_density.csv"),
        by_domain_ub = joinpath(output_dir, "exact_res_by_domain_ub.csv"),
        by_max_distinct_coeffs = joinpath(output_dir, "exact_res_by_max_distinct_coeffs.csv"),
    )
end

function write_aggregate_csvs(output_dir::AbstractString, aggregate_rows)
    mkpath(output_dir)
    paths = aggregate_output_paths(output_dir)
    write_csv(paths.by_nvars, aggregate_rows.by_nvars)
    write_csv(paths.by_density, aggregate_rows.by_density)
    write_csv(paths.by_domain_ub, aggregate_rows.by_domain_ub)
    write_csv(paths.by_max_distinct_coeffs, aggregate_rows.by_max_distinct_coeffs)
    return paths
end

function print_config(result)
    config = result.config
    println("Exact bound reduction family experiment")
    println("count = $(config.count)")
    println("nvars = $(join(config.nvars, ","))")
    println("domain_lb = $(config.domain_lb)")
    println("domain_ubs = $(join(config.domain_ubs, ","))")
    println("max_distinct_coeffs = $(join(config.max_distinct_coeffs, ","))")
    println("seed_base = $(config.seed_base)")
    println("seed_step = $(config.seed_step)")
    println("densities = $(join(config.densities, ","))")
    println("coeff_range = $(config.coeff_lb):$(config.coeff_ub) excluding 0")
    println("offset_range = $(config.offset_lb):$(config.offset_ub)")
    config.output_dir !== nothing && println("output_dir = $(config.output_dir)")
    println("parameter_combinations = $(length(result.results))")
    println()
    return nothing
end

function _format_group_value(value)
    return value isa AbstractFloat ? @sprintf("%.6f", value) : string(value)
end

function print_aggregate_table(title::AbstractString, rows, field::Symbol)
    println(title)
    @printf(
        "%20s %12s %14s %14s %10s %10s %10s %10s %14s\n",
        string(field),
        "constraints",
        "opt_avg_red",
        "opt_std_red",
        "tw_1",
        "tw_2",
        "tw_3",
        "tw_ge4",
        "avg_time_sec",
    )
    for row in rows
        @printf(
            "%20s %12d %14.6f %14.6f %10.2f %10.2f %10.2f %10.2f %14.6f\n",
            _format_group_value(getproperty(row, field)),
            row.constraints,
            row.opt_avg_red,
            row.opt_std_red,
            row.tw_1,
            row.tw_2,
            row.tw_3,
            row.tw_ge4,
            row.avg_wall_time_sec_per_constraint,
        )
    end
    return nothing
end

function print_aggregate_tables(aggregate_rows)
    print_aggregate_table("Grouped by nvars", aggregate_rows.by_nvars, :nvars)
    println()
    print_aggregate_table("Grouped by density", aggregate_rows.by_density, :density)
    println()
    print_aggregate_table("Grouped by domain_ub", aggregate_rows.by_domain_ub, :domain_ub)
    println()
    print_aggregate_table(
        "Grouped by max_distinct_coeffs",
        aggregate_rows.by_max_distinct_coeffs,
        :max_distinct_coeffs,
    )
    return nothing
end

function print_written_paths(paths)
    println()
    println("Wrote CSVs:")
    for name in propertynames(paths)
        println("  $(getproperty(paths, name))")
    end
    return nothing
end

function main(args::Vector{String} = copy(ARGS))
    config = build_config(args)
    config === nothing && return nothing

    result = run_experiment(config)
    print_config(result)
    print_aggregate_tables(result.aggregate_rows)

    output_paths = nothing
    if config.output_dir !== nothing
        output_paths = write_aggregate_csvs(config.output_dir, result.aggregate_rows)
        print_written_paths(output_paths)
    end

    return merge(result, (output_paths = output_paths,))
end

end # module

if _EXACT_BOUND_REDUCTION_FAMILY_RUNNING_AS_SCRIPT
    ExactBoundReductionFamilyExperiment.main(copy(ARGS))
end
