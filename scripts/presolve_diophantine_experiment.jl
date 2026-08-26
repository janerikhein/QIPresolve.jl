#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
end

module DiophantineParityPresolveExperiment

using JuMP: backend
using Printf
using Random

import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC
using QIPresolve.InstanceGeneration: generate_random_qip_model

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}

const DEFAULT_N_VARS = 100
const DEFAULT_INSTANCES_PER_CELL = 100
const DEFAULT_CONSTRAINT_COUNTS = [50, 75, 100, 125, 150, 175]
const DEFAULT_TYPES = ["bilinear", "separable", "pure", "general"]
const DEFAULT_SEED_BASE = 1
const DEFAULT_SEED_STEP = 1
const DEFAULT_OUTPUT_DIR = joinpath("results", "diophantine_presolve_experiment")
const DEFAULT_PARITY_STRATEGY = QIP.PresolveConfig.DEFAULT_PRESOLVE_PARITY_STRATEGY
const SUMMARY_FILENAME = "parity_presolve_summary.csv"
const SUMMARY_HEADER = "type,m,avg_b_fix,avg_b_pat,avg_dom_red,avg_time"

const CLI_KEYS = Dict(
    "nvars" => :nvars,
    "n-vars" => :nvars,
    "instances-per-cell" => :instances_per_cell,
    "instances" => :instances_per_cell,
    "count" => :instances_per_cell,
    "constraints" => :constraints,
    "constraint-counts" => :constraints,
    "m" => :constraints,
    "types" => :types,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "output-dir" => :output_dir,
    "output" => :output_dir,
    "parity-strategy" => :parity_strategy,
)

Base.@kwdef struct CliConfig
    nvars::Int = DEFAULT_N_VARS
    instances_per_cell::Int = DEFAULT_INSTANCES_PER_CELL
    constraints::Vector{Int} = copy(DEFAULT_CONSTRAINT_COUNTS)
    types::Vector{String} = copy(DEFAULT_TYPES)
    seed_base::Int = DEFAULT_SEED_BASE
    seed_step::Int = DEFAULT_SEED_STEP
    output_dir::String = abspath(DEFAULT_OUTPUT_DIR)
    parity_strategy::Symbol = DEFAULT_PARITY_STRATEGY
end

function usage()
    return """
    Usage:
      julia --project=. scripts/presolve_diophantine_experiment.jl [options]

    Options:
      --nvars n                    Original variables per instance, default $DEFAULT_N_VARS
      --instances-per-cell n       Instances per (type, m), default $DEFAULT_INSTANCES_PER_CELL
      --constraints list           Comma-separated constraint counts, default $(join(DEFAULT_CONSTRAINT_COUNTS, ","))
      --types list                 Comma-separated types, default $(join(DEFAULT_TYPES, ","))
      --seed-base n                First deterministic instance seed, default $DEFAULT_SEED_BASE
      --seed-step n                Seed increment, default $DEFAULT_SEED_STEP
      --output-dir path            Output directory, default $DEFAULT_OUTPUT_DIR
      --parity-strategy name       full, mod2-basic, or mod4-basic, default $DEFAULT_PARITY_STRATEGY
      -h, --help                   Show this help
    """
end

function parse_int(value::AbstractString, name::AbstractString)::Int
    parsed = tryparse(Int, strip(value))
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function parse_int_list(value::AbstractString, name::AbstractString)::Vector{Int}
    values = Int[]
    for part in split(value, ","; keepempty = false)
        push!(values, parse_int(part, name))
    end
    isempty(values) && error("$name must contain at least one integer")
    return values
end

function canonical_type(value::AbstractString)::String
    normalized = lowercase(strip(value))
    normalized = replace(normalized, "_" => "-")
    normalized == "seperable" && return "separable"
    normalized == "separable-quadratic" && return "separable"
    normalized == "purely-quadratic" && return "pure"
    normalized == "pure-quadratic" && return "pure"
    normalized in DEFAULT_TYPES && return normalized
    error("Invalid type: $value. Expected one of $(join(DEFAULT_TYPES, ", ")).")
end

function parse_types(value::AbstractString)::Vector{String}
    types = String[]
    seen = Set{String}()
    for part in split(value, ","; keepempty = false)
        type = canonical_type(part)
        type in seen && continue
        push!(types, type)
        push!(seen, type)
    end
    isempty(types) && error("types must contain at least one type")
    return types
end

function parse_parity_strategy(value::AbstractString)::Symbol
    normalized = Symbol(replace(lowercase(strip(value)), "-" => "_"))
    normalized in (:full, :mod2_basic, :mod4_basic) && return normalized
    error("Invalid parity strategy: $value. Expected full, mod2-basic, or mod4-basic.")
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

function build_config(args::Vector{String})::Union{Nothing, CliConfig}
    options = parse_raw_options(args)
    options === nothing && return nothing

    config = CliConfig(
        nvars = parse_int(get(options, :nvars, string(DEFAULT_N_VARS)), "nvars"),
        instances_per_cell = parse_int(
            get(options, :instances_per_cell, string(DEFAULT_INSTANCES_PER_CELL)),
            "instances_per_cell",
        ),
        constraints = haskey(options, :constraints) ?
            parse_int_list(options[:constraints], "constraints") :
            copy(DEFAULT_CONSTRAINT_COUNTS),
        types = haskey(options, :types) ? parse_types(options[:types]) : copy(DEFAULT_TYPES),
        seed_base = parse_int(get(options, :seed_base, string(DEFAULT_SEED_BASE)), "seed_base"),
        seed_step = parse_int(get(options, :seed_step, string(DEFAULT_SEED_STEP)), "seed_step"),
        output_dir = abspath(get(options, :output_dir, DEFAULT_OUTPUT_DIR)),
        parity_strategy = haskey(options, :parity_strategy) ?
            parse_parity_strategy(options[:parity_strategy]) :
            DEFAULT_PARITY_STRATEGY,
    )

    config.nvars >= 1 || error("nvars must be >= 1")
    config.instances_per_cell >= 1 || error("instances_per_cell must be >= 1")
    all(>=(0), config.constraints) || error("constraint counts must be nonnegative")
    config.seed_base >= 1 || error("seed_base must be >= 1")
    config.seed_step >= 0 || error("seed_step must be >= 0")
    if "bilinear" in config.types
        config.nvars >= 2 || error("bilinear instances require nvars >= 2")
    end

    return config
end

function type_probabilities(type::AbstractString, p::Real)
    canonical = canonical_type(type)
    density = Float64(p)
    0.0 <= density <= 1.0 || error("p must be in [0, 1], got $p")

    if canonical == "bilinear"
        return (p_var_bilin = density, p_var_diag = 0.0, p_var_lin = 0.0)
    elseif canonical == "separable"
        return (p_var_bilin = 0.0, p_var_diag = density, p_var_lin = density)
    elseif canonical == "pure"
        return (p_var_bilin = density, p_var_diag = density, p_var_lin = 0.0)
    elseif canonical == "general"
        return (p_var_bilin = density, p_var_diag = density, p_var_lin = density)
    end

    error("Invalid type: $type")
end

function generator_kwargs(type::AbstractString, p::Real)
    probs = type_probabilities(type, p)
    return (
        p_con_eq = 1.0,
        var_threshold_lb = 0,
        var_threshold_ub = 100,
        p_var_is_candidate = 1.0,
        p_var_bilin = probs.p_var_bilin*probs.p_var_bilin,
        p_var_diag = probs.p_var_diag,
        p_var_lin = probs.p_var_lin,
        coeff_lb = -10,
        coeff_ub = 10,
        force_diag_even = false,
        force_lin_even = false,
        force_feasibility = true,
        constraint_slack_range = [0],
    )
end

function _empty_objective()
    return PC.QuadExpr(QuadTerm[], LinTerm[])
end

function _sample_density(rng::AbstractRNG)::Float64

    while true
        p = 0.19 * rand(rng) + 0.01
        0.01 < p < 0.2 && return p
    end
end

function _sample_bounds_and_reference(rng::AbstractRNG, nvars::Int)
    threshold = rand(rng, 5:100)
    upper_bounds = [rand(rng, 1:threshold) for _ in 1:nvars]
    x_star = [rand(rng, 0:upper_bounds[var_id]) for var_id in 1:nvars]
    return upper_bounds, x_star
end

function postprocess_diophantine_system!(
        model::PC.QPModel,
        upper_bounds::AbstractVector{<:Integer},
        x_star::AbstractVector{<:Integer},
    )
    model.obj_expr = _empty_objective()
    model.obj_sense = :min

    for var_id in eachindex(upper_bounds)
        PC.set_var_bounds!(model, var_id, 0.0, Float64(upper_bounds[var_id]))
    end
    model._max_var_id = max(model._max_var_id, length(upper_bounds))

    for con in model.cons
        rhs = PC.eval_full(con.qe, x_star)
        con.lhs = rhs
        con.rhs = rhs
        PC.normalize!(con)
    end

    return model
end

function build_random_diophantine_model(
        type::AbstractString,
        nvars::Int,
        ncons::Int,
        p::Real,
        generator_seed::Int,
        postprocess_rng::AbstractRNG,
    )
    jump_model, _ = generate_random_qip_model(
        nvars,
        ncons;
        generator_kwargs(type, p)...,
        seed = generator_seed,
    )
    model = QIP.build_model(QIP.from_moi(backend(jump_model)))
    upper_bounds, x_star = _sample_bounds_and_reference(postprocess_rng, nvars)
    postprocess_diophantine_system!(model, upper_bounds, x_star)
    return model, upper_bounds, x_star
end

log_domain_sum(model::PC.QPModel) = sum(log(var.ub - var.lb + 1.0) for var in values(model.vars); init = 0.0)

function average_domain_size(upper_bounds::AbstractVector{<:Integer})::Float64
    isempty(upper_bounds) && return 0.0
    return sum(upper + 1 for upper in upper_bounds; init = 0) / length(upper_bounds)
end

function total_original_bits(upper_bounds::AbstractVector{<:Integer})::Int
    return sum(ndigits(upper; base = 2) for upper in upper_bounds; init = 0)
end

function _active_var_is_fixed(model::PC.QPModel, var_id::PC.VarId)::Bool
    haskey(model.vars, var_id) || return false
    var = model.vars[var_id]
    return var.lb == var.ub
end

function _stored_bit_is_fixed(
        bit::PC.StoredBit,
        postsolver::PC.ParityPostsolver,
        model::PC.QPModel,
        cache::Dict{PC.VarId, Bool},
        active::Set{PC.VarId},
    )::Bool
    (bit.kind == PC.FIXED0 || bit.kind == PC.FIXED1) && return true

    ref_var_id = bit.ref_var_id
    ref_var_id === nothing && return false
    return _reconstruction_is_fixed(postsolver, model, ref_var_id, cache, active)
end

function _high_order_is_fixed(
        postsolver::PC.ParityPostsolver,
        model::PC.QPModel,
        owner_var_id::PC.VarId,
        var_data::PC.VarReconstruction,
        cache::Dict{PC.VarId, Bool},
        active::Set{PC.VarId},
    )::Bool
    current_var_id = var_data.current_var_id
    if current_var_id === nothing
        return var_data.fixed_high_order !== nothing
    elseif current_var_id == owner_var_id
        return _active_var_is_fixed(model, current_var_id)
    elseif haskey(postsolver.var_data, current_var_id)
        return _reconstruction_is_fixed(postsolver, model, current_var_id, cache, active)
    else
        return _active_var_is_fixed(model, current_var_id)
    end
end

function _reconstruction_is_fixed(
        postsolver::PC.ParityPostsolver,
        model::PC.QPModel,
        var_id::PC.VarId,
        cache::Dict{PC.VarId, Bool},
        active::Set{PC.VarId},
    )::Bool
    haskey(cache, var_id) && return cache[var_id]
    if !(haskey(postsolver.var_data, var_id))
        fixed = _active_var_is_fixed(model, var_id)
        cache[var_id] = fixed
        return fixed
    end
    var_id in active && return false

    push!(active, var_id)
    var_data = postsolver.var_data[var_id]
    fixed = _high_order_is_fixed(postsolver, model, var_id, var_data, cache, active)
    if fixed
        for bit in var_data.bits
            if !_stored_bit_is_fixed(bit, postsolver, model, cache, active)
                fixed = false
                break
            end
        end
    end
    delete!(active, var_id)
    cache[var_id] = fixed
    return fixed
end

function postsolve_bit_counts(
        postsolver::PC.ParityPostsolver,
        model::PC.QPModel,
        upper_bounds::AbstractVector{<:Integer},
    )
    fixed_bits = 0
    pattern_bits = 0
    fixed_cache = Dict{PC.VarId, Bool}()

    for var_id in 1:length(upper_bounds)
        width = ndigits(upper_bounds[var_id]; base = 2)
        var_data = get(postsolver.var_data, var_id, nothing)
        if var_data === nothing
            continue
        end

        stored_count = min(width, length(var_data.bits))
        for bit_index in 1:stored_count
            bit = var_data.bits[bit_index]
            if _stored_bit_is_fixed(bit, postsolver, model, fixed_cache, Set{PC.VarId}())
                fixed_bits += 1
            elseif bit.kind == PC.BINVAR || bit.kind == PC.NEGATED_BINVAR
                pattern_bits += 1
            end
        end

        remaining_bits = max(0, width - length(var_data.bits))
        if remaining_bits > 0 &&
                _high_order_is_fixed(postsolver, model, var_id, var_data, fixed_cache, Set{PC.VarId}())
            fixed_bits += remaining_bits
        end
    end

    total_bits = total_original_bits(upper_bounds)
    return (
        fixed_bits = fixed_bits,
        pattern_bits = pattern_bits,
        total_bits = total_bits,
        b_fix = total_bits == 0 ? 0.0 : fixed_bits / total_bits,
        b_pat = total_bits == 0 ? 0.0 : pattern_bits / total_bits,
    )
end

function run_instance(
        type::AbstractString,
        nvars::Int,
        ncons::Int,
        seed::Int;
        parity_strategy::Symbol = DEFAULT_PARITY_STRATEGY,
    )
    rng = MersenneTwister(seed)
    p = _sample_density(rng)
    generator_seed = rand(rng, 1:typemax(Int32))
    model, upper_bounds, _ = build_random_diophantine_model(
        type,
        nvars,
        ncons,
        p,
        generator_seed,
        rng,
    )

    original_log_domain_sum = log_domain_sum(model)
    postsolver = PC.ParityPostsolver(keys(model.vars))
    propagator = PC.PropagationManager(PC.VarId[])
    PC.normalize!(model, postsolver)
    @assert model.infeasible==false
    start_time = time()
    PC.parity_presolve!(
        model,
        propagator,
        postsolver,
        nothing;
        parity_strategy = parity_strategy,
    )
    presolve_time = time() - start_time
    @assert model.infeasible==false
    presolved_log_domain_sum = model.infeasible ? 0.0 : log_domain_sum(model)
    dom_red = original_log_domain_sum == 0.0 ?
        0.0 :
        (original_log_domain_sum - presolved_log_domain_sum) / original_log_domain_sum
    bit_counts = postsolve_bit_counts(postsolver, model, upper_bounds)

    return (
        type = canonical_type(type),
        ncons = ncons,
        seed = seed,
        p = p,
        avg_domain_size = average_domain_size(upper_bounds),
        b_fix = bit_counts.b_fix,
        b_pat = bit_counts.b_pat,
        dom_red = dom_red,
        time = presolve_time,
        infeasible = model.infeasible,
    )
end

function summary_path(config::CliConfig)
    return joinpath(config.output_dir, SUMMARY_FILENAME)
end

function density_path(config::CliConfig, type::AbstractString)
    return joinpath(config.output_dir, "density_vs_dom_red_$(canonical_type(type)).txt")
end

function avg_domain_path(config::CliConfig, type::AbstractString)
    return joinpath(config.output_dir, "avg_domain_size_vs_dom_red_$(canonical_type(type)).txt")
end

function prepare_output_files!(config::CliConfig)
    mkpath(config.output_dir)
    open(summary_path(config), "w") do io
        println(io, SUMMARY_HEADER)
    end

    for type in config.types
        open(density_path(config, type), "w") do _ end
        open(avg_domain_path(config, type), "w") do _ end
    end

    return nothing
end

function format_float(value::Real)::String
    return @sprintf("%.17g", Float64(value))
end

function format_summary_row(row)::String
    return join(
        (
            row.type,
            string(row.m),
            format_float(row.avg_b_fix),
            format_float(row.avg_b_pat),
            format_float(row.avg_dom_red),
            format_float(row.avg_time),
        ),
        ",",
    )
end

function append_summary_row!(config::CliConfig, row)
    open(summary_path(config), "a") do io
        println(io, format_summary_row(row))
    end
    return nothing
end

function append_tuple_row!(path::AbstractString, x::Real, y::Real)
    open(path, "a") do io
        @printf(io, "(%.17g, %.17g)\n", Float64(x), Float64(y))
    end
    return nothing
end

function cell_summary(type::String, ncons::Int, instances)
    count = length(instances)
    count > 0 || error("Cannot summarize an empty cell")
    return (
        type = type,
        m = ncons,
        avg_b_fix = sum(instance.b_fix for instance in instances; init = 0.0) / count,
        avg_b_pat = sum(instance.b_pat for instance in instances; init = 0.0) / count,
        avg_dom_red = sum(instance.dom_red for instance in instances; init = 0.0) / count,
        avg_time = sum(instance.time for instance in instances; init = 0.0) / count,
    )
end

function run_experiment(config::CliConfig)
    prepare_output_files!(config)
    rows = NamedTuple[]
    instance_counter = 0

    for type in config.types
        for ncons in config.constraints
            instances = NamedTuple[]
            sizehint!(instances, config.instances_per_cell)

            for _ in 1:config.instances_per_cell
                seed = config.seed_base + instance_counter * config.seed_step
                instance_counter += 1
                instance = run_instance(
                    type,
                    config.nvars,
                    ncons,
                    seed;
                    parity_strategy = config.parity_strategy,
                )
                push!(instances, instance)
                append_tuple_row!(density_path(config, type), instance.p, instance.dom_red)
                append_tuple_row!(
                    avg_domain_path(config, type),
                    instance.avg_domain_size,
                    instance.dom_red,
                )
            end

            row = cell_summary(type, ncons, instances)
            push!(rows, row)
            append_summary_row!(config, row)
            println(format_summary_row(row))
            flush(stdout)
        end
    end

    return (
        config = config,
        rows = rows,
        summary_path = summary_path(config),
        instance_count = instance_counter,
    )
end

function main(args::Vector{String} = copy(ARGS))
    config = build_config(args)
    config === nothing && return nothing
    return run_experiment(config)
end

end # module

if _RUNNING_AS_SCRIPT
    DiophantineParityPresolveExperiment.main(copy(ARGS))
end
