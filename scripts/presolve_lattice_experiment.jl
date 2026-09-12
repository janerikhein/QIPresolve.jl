#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
end

module LatticeParityPresolveExperiment

using JuMP: backend
using Printf

import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC
using QIPresolve.InstanceGeneration:
    generate_2_connected_instance,
    generate_globally_rigid_instance,
    generate_laman_instance

const DEFAULT_FREE_VERTICES = 50
const DEFAULT_R = 50
const DEFAULT_NUM_ANCHORS = 0
const DEFAULT_INSTANCES_PER_TYPE = 300
const DEFAULT_TYPES = [
    "2-connected-sparse",
    "2-connected-dense",
    "laman",
    "globally-rigid",
]
const DEFAULT_SEED_BASE = 1
const DEFAULT_SEED_STEP = 1
const DEFAULT_OUTPUT_DIR = joinpath("results", "lattice_presolve_experiment")
const DEFAULT_PARITY_STRATEGY = QIP.PresolveConfig.DEFAULT_PRESOLVE_PARITY_STRATEGY
const SUMMARY_FILENAME = "lattice_parity_presolve_summary.csv"
const SUMMARY_HEADER = "type,m,avg_b_fix,avg_b_pat,avg_dom_red,avg_time"

const CLI_KEYS = Dict(
    "num-anchors" => :num_anchors,
    "anchors" => :num_anchors,
    "instances-per-type" => :instances_per_type,
    "instances" => :instances_per_type,
    "count" => :instances_per_type,
    "types" => :types,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "output-dir" => :output_dir,
    "output" => :output_dir,
    "parity-strategy" => :parity_strategy,
    "free-vertices" => :free_vertices,
)

Base.@kwdef struct CliConfig
    num_anchors::Int = DEFAULT_NUM_ANCHORS
    instances_per_type::Int = DEFAULT_INSTANCES_PER_TYPE
    types::Vector{String} = copy(DEFAULT_TYPES)
    seed_base::Int = DEFAULT_SEED_BASE
    seed_step::Int = DEFAULT_SEED_STEP
    output_dir::String = abspath(DEFAULT_OUTPUT_DIR)
    parity_strategy::Symbol = DEFAULT_PARITY_STRATEGY
    free_vertices::Int = DEFAULT_FREE_VERTICES
end

function usage()
    return """
    Usage:
      julia --project=. scripts/presolve_lattice_experiment.jl [options]

    Options:
      --num-anchors n             Number of anchored vertices, default $DEFAULT_NUM_ANCHORS
      --instances-per-type n      Instances per lattice type, default $DEFAULT_INSTANCES_PER_TYPE
      --types list                Comma-separated types, default $(join(DEFAULT_TYPES, ","))
      --seed-base n               First deterministic instance seed, default $DEFAULT_SEED_BASE
      --seed-step n               Seed increment, default $DEFAULT_SEED_STEP
      --output-dir path           Output directory, default $DEFAULT_OUTPUT_DIR
      --parity-strategy name      full, mod2-basic, or mod4-basic, default $DEFAULT_PARITY_STRATEGY
      --free-vertices n           Free vertices used to derive n, default $DEFAULT_FREE_VERTICES
      -h, --help                  Show this help
    """
end

function parse_int(value::AbstractString, name::AbstractString)::Int
    parsed = tryparse(Int, strip(value))
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function canonical_type(value::AbstractString)::String
    normalized = lowercase(strip(value))
    normalized = replace(normalized, "_" => "-", " " => "-")

    if normalized in ("2-connected-sparse", "two-connected-sparse", "con-sparse", "sparse")
        return "2-connected-sparse"
    elseif normalized in ("2-connected-dense", "two-connected-dense", "con-dense", "dense")
        return "2-connected-dense"
    elseif normalized in ("laman", "lam")
        return "laman"
    elseif normalized in ("globally-rigid", "global-rigid", "globally-rigid-graph", "gr")
        return "globally-rigid"
    end

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

effective_n(num_anchors::Int, free_vertices::Int = DEFAULT_FREE_VERTICES)::Int =
    num_anchors == 0 ? free_vertices + 1 : free_vertices + num_anchors

effective_n(config::CliConfig)::Int = effective_n(config.num_anchors, config.free_vertices)

function _minimum_n(type::AbstractString)::Int
    canonical = canonical_type(type)
    canonical == "globally-rigid" && return 4
    return 3
end

function build_config(args::Vector{String})::Union{Nothing, CliConfig}
    options = parse_raw_options(args)
    options === nothing && return nothing

    config = CliConfig(
        num_anchors = parse_int(get(options, :num_anchors, string(DEFAULT_NUM_ANCHORS)), "num_anchors"),
        instances_per_type = parse_int(
            get(options, :instances_per_type, string(DEFAULT_INSTANCES_PER_TYPE)),
            "instances_per_type",
        ),
        types = haskey(options, :types) ? parse_types(options[:types]) : copy(DEFAULT_TYPES),
        seed_base = parse_int(get(options, :seed_base, string(DEFAULT_SEED_BASE)), "seed_base"),
        seed_step = parse_int(get(options, :seed_step, string(DEFAULT_SEED_STEP)), "seed_step"),
        output_dir = abspath(get(options, :output_dir, DEFAULT_OUTPUT_DIR)),
        parity_strategy = haskey(options, :parity_strategy) ?
            parse_parity_strategy(options[:parity_strategy]) :
            DEFAULT_PARITY_STRATEGY,
        free_vertices = parse_int(
            get(options, :free_vertices, string(DEFAULT_FREE_VERTICES)),
            "free_vertices",
        ),
    )

    config.num_anchors >= 0 || error("num_anchors must be >= 0")
    config.instances_per_type >= 1 || error("instances_per_type must be >= 1")
    config.seed_base >= 1 || error("seed_base must be >= 1")
    config.seed_step >= 0 || error("seed_step must be >= 0")
    config.free_vertices >= 0 || error("free_vertices must be >= 0")

    n = effective_n(config)
    for type in config.types
        n >= _minimum_n(type) || error(
            "$type instances require effective n >= $(_minimum_n(type)), got $n"
        )
    end

    return config
end

function _edge_density(type::AbstractString)::Float64
    canonical = canonical_type(type)
    canonical == "2-connected-sparse" && return 0.05
    canonical == "2-connected-dense" && return 0.3
    error("Lattice type $canonical does not use edge density")
end

function lattice_edge_count(type::AbstractString, n::Int)::Int
    canonical = canonical_type(type)
    if canonical in ("2-connected-sparse", "2-connected-dense")
        m_min = n
        m_max = n * (n - 1) ÷ 2
        return round(Int, m_min + _edge_density(canonical) * (m_max - m_min))
    elseif canonical == "laman"
        return 2n - 3
    elseif canonical == "globally-rigid"
        return 2n - 2
    end

    error("Invalid type: $type")
end

function build_lattice_jump_model(
        type::AbstractString,
        n::Int,
        seed::Int,
        num_anchors::Int,
    )
    canonical = canonical_type(type)
    if canonical in ("2-connected-sparse", "2-connected-dense")
        return generate_2_connected_instance(
            n;
            R = DEFAULT_R,
            edge_density = _edge_density(canonical),
            seed = seed,
            num_anchors = num_anchors,
            alpha = 0.0,
        )
    elseif canonical == "laman"
        return generate_laman_instance(
            n;
            R = DEFAULT_R,
            seed = seed,
            num_anchors = num_anchors,
            alpha = 0.0,
        )
    elseif canonical == "globally-rigid"
        return generate_globally_rigid_instance(
            n;
            R = DEFAULT_R,
            seed = seed,
            num_anchors = num_anchors,
            alpha = 0.0,
        )
    end

    error("Invalid type: $type")
end

function build_lattice_model(
        type::AbstractString,
        n::Int,
        seed::Int,
        num_anchors::Int,
    )::PC.QPModel
    jump_model, _, _ = build_lattice_jump_model(type, n, seed, num_anchors)
    return QIP.build_model(QIP.from_moi(backend(jump_model)))
end

log_domain_sum(model::PC.QPModel) =
    sum(log(var.ub - var.lb + 1.0) for var in values(model.vars); init = 0.0)

function bit_width(var::PC.IntVar)::Int
    isfinite(var.lb) && isfinite(var.ub) || error(
        "Cannot compute bit width for non-finite domain [$(var.lb), $(var.ub)]"
    )
    lower = ceil(Int, var.lb)
    upper = floor(Int, var.ub)
    lower <= upper || error("Cannot compute bit width for empty domain [$lower, $upper]")

    return max(1, ndigits(upper - lower; base = 2))
end

function original_bit_widths(model::PC.QPModel)::Dict{PC.VarId, Int}
    return Dict(var_id => bit_width(var) for (var_id, var) in model.vars)
end

total_original_bits(bit_widths::AbstractDict{PC.VarId, Int}) =
    sum(values(bit_widths); init = 0)

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
    if !haskey(postsolver.var_data, var_id)
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
        bit_widths::AbstractDict{PC.VarId, Int},
    )
    fixed_bits = 0
    pattern_bits = 0
    fixed_cache = Dict{PC.VarId, Bool}()

    for var_id in sort!(collect(keys(bit_widths)))
        width = bit_widths[var_id]
        var_data = get(postsolver.var_data, var_id, nothing)
        var_data === nothing && continue

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

    total_bits = total_original_bits(bit_widths)
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
        n::Int,
        seed::Int,
        num_anchors::Int;
        parity_strategy::Symbol = DEFAULT_PARITY_STRATEGY,
    )
    canonical = canonical_type(type)
    model = build_lattice_model(canonical, n, seed, num_anchors)

    original_log_domain_sum = log_domain_sum(model)
    bit_widths = original_bit_widths(model)
    postsolver = PC.ParityPostsolver(keys(model.vars))
    propagator = PC.PropagationManager(PC.VarId[])

    PC.normalize!(model, postsolver)
    model.infeasible && error("Generated $canonical instance with seed $seed became infeasible before parity presolve")

    start_time = time()
    PC.parity_presolve!(
        model,
        propagator,
        postsolver,
        nothing;
        parity_strategy = parity_strategy,
    )
    presolve_time = time() - start_time
    model.infeasible && error("Generated $canonical instance with seed $seed became infeasible during parity presolve")

    presolved_log_domain_sum = log_domain_sum(model)
    dom_red = original_log_domain_sum == 0.0 ?
        0.0 :
        (original_log_domain_sum - presolved_log_domain_sum) / original_log_domain_sum
    bit_counts = postsolve_bit_counts(postsolver, model, bit_widths)

    return (
        type = canonical,
        m = lattice_edge_count(canonical, n),
        n = n,
        seed = seed,
        num_anchors = num_anchors,
        b_fix = bit_counts.b_fix,
        b_pat = bit_counts.b_pat,
        dom_red = dom_red,
        time = presolve_time,
    )
end

summary_path(config::CliConfig) = joinpath(config.output_dir, SUMMARY_FILENAME)

function prepare_output_files!(config::CliConfig)
    mkpath(config.output_dir)
    open(summary_path(config), "w") do io
        println(io, SUMMARY_HEADER)
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

function cell_summary(type::String, m::Int, instances)
    count = length(instances)
    count > 0 || error("Cannot summarize an empty cell")
    return (
        type = type,
        m = m,
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
    n = effective_n(config)

    for type in config.types
        instances = NamedTuple[]
        sizehint!(instances, config.instances_per_type)

        for _ in 1:config.instances_per_type
            seed = config.seed_base + instance_counter * config.seed_step
            instance_counter += 1
            instance = run_instance(
                type,
                n,
                seed,
                config.num_anchors;
                parity_strategy = config.parity_strategy,
            )
            push!(instances, instance)
        end

        row = cell_summary(type, lattice_edge_count(type, n), instances)
        push!(rows, row)
        append_summary_row!(config, row)
        println(format_summary_row(row))
        flush(stdout)
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
    LatticeParityPresolveExperiment.main(copy(ARGS))
end
