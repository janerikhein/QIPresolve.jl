#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module MergeGraphEmbeddingPresolveStats

import CSV
import JSON

const DEFAULT_BASE_DIR = joinpath("results", "graph_embedding_instances")
const DEFAULT_INSTANCES_CSV = joinpath(DEFAULT_BASE_DIR, "instances.csv")
const DEFAULT_STATS_DIR = joinpath(DEFAULT_BASE_DIR, "presolve_stats")
const DEFAULT_OUTPUT = joinpath(DEFAULT_BASE_DIR, "instances_with_presolve_stats.csv")

const INSTANCE_COLUMNS = [
    "instance_name",
    "type",
    "n",
    "R",
    "alpha",
    "edge_density",
    "infeas_strategy",
    "infeas_base",
    "box_scale",
]

const STATS_COLUMNS = [
    "log_domain_sum_orig",
    "log_domain_sum_ps",
    "num_quadratic_constraints_orig",
    "num_quadratic_constraints_ps",
    "num_linear_constraints_orig",
    "num_linear_constraints_ps",
    "parity_presolve_time",
    "residue_presolve_time",
    "num_parity_substitutions",
    "num_pattern_substitutions",
    "num_variables_fixed",
    "num_constraints_tightened",
    "num_constraints_eq_tightened",
    "avg_relative_bound_tightening",
    "infeasibility_detected",
    "infeasibility_source",
]

const OUTPUT_COLUMNS = [INSTANCE_COLUMNS; STATS_COLUMNS]

Base.@kwdef struct CliConfig
    instances_csv::String = abspath(DEFAULT_INSTANCES_CSV)
    stats_dir::String = abspath(DEFAULT_STATS_DIR)
    output::String = abspath(DEFAULT_OUTPUT)
end

function usage()
    return """
    Usage:
      julia --project=. scripts/merge_graph_embedding_presolve_stats.jl [options]

    Options:
      --instances-csv path    Input graph embedding instances CSV
      --stats-dir dir         Directory containing <instance_name>_stats.json files
      --output path           Output merged CSV path
      -h, --help              Show this help
    """
end

function _option_value(args::Vector{String}, index::Int, option::String)
    if occursin("=", option)
        parts = split(option, "="; limit = 2)
        return parts[1], parts[2], 1
    end

    index < length(args) || error("Missing value for option $option")
    return option, args[index + 1], 2
end

function parse_args(args::Vector{String})
    if !isempty(args) && args[1] in ("-h", "--help")
        println(usage())
        return nothing
    end

    instances_csv = DEFAULT_INSTANCES_CSV
    stats_dir = DEFAULT_STATS_DIR
    output = DEFAULT_OUTPUT

    index = 1
    while index <= length(args)
        arg = args[index]
        startswith(arg, "--") || error("Unexpected positional argument: $arg")

        option, value, consumed = _option_value(args, index, arg)
        if option == "--instances-csv"
            instances_csv = value
        elseif option == "--stats-dir"
            stats_dir = value
        elseif option == "--output"
            output = value
        else
            error("Unknown option: $option")
        end
        index += consumed
    end

    return CliConfig(
        instances_csv = abspath(instances_csv),
        stats_dir = abspath(stats_dir),
        output = abspath(output),
    )
end

function _require_key(dict::AbstractDict, context::AbstractString, key::AbstractString)
    haskey(dict, key) || error("Missing key $context.$key")
    return dict[key]
end

function _require_section(data::AbstractDict, section::AbstractString)
    value = _require_key(data, "stats JSON", section)
    value isa AbstractDict || error("Expected stats JSON.$section to be an object")
    return value
end

function _row_value(row, column::AbstractString)
    value = getproperty(row, Symbol(column))
    return value === missing ? missing : value
end

function _ordered_named_tuple(values::AbstractDict, columns::Vector{String})
    names = Tuple(Symbol.(columns))
    return NamedTuple{names}(Tuple(values[column] for column in columns))
end

function _stats_json_path(stats_dir::AbstractString, instance_name::AbstractString)
    return joinpath(stats_dir, "$(instance_name)_stats.json")
end

function _infeasibility_source(parity_stats::AbstractDict, residue_stats::AbstractDict)
    parity_source = string(_require_key(parity_stats, "parity_stats", "infeasibility_source"))
    parity_source != "none" && return parity_source

    residue_detected = Bool(_require_key(residue_stats, "residue_stats", "infeasibility_detected"))
    return residue_detected ? "residue" : "none"
end

function _stats_values(data::AbstractDict, instance_name::AbstractString)
    json_instance_name = string(_require_key(data, "stats JSON", "instance_name"))
    json_instance_name == instance_name || error(
        "Stats JSON instance_name mismatch for $instance_name: found $json_instance_name"
    )

    model_loaded = _require_section(data, "model_loaded")
    model_after = _require_section(data, "model_after")
    parity_stats = _require_section(data, "parity_stats")
    residue_stats = _require_section(data, "residue_stats")

    parity_infeasible = Bool(_require_key(parity_stats, "parity_stats", "infeasibility_detected"))
    residue_infeasible = Bool(_require_key(residue_stats, "residue_stats", "infeasibility_detected"))

    return Dict{String, Any}(
        "log_domain_sum_orig" => _require_key(model_loaded, "model_loaded", "log_domain_sum"),
        "log_domain_sum_ps" => _require_key(model_after, "model_after", "log_domain_sum"),
        "num_quadratic_constraints_orig" =>
            _require_key(model_loaded, "model_loaded", "num_quadratic_constraints"),
        "num_quadratic_constraints_ps" =>
            _require_key(model_after, "model_after", "num_quadratic_constraints"),
        "num_linear_constraints_orig" =>
            _require_key(model_loaded, "model_loaded", "num_linear_constraints"),
        "num_linear_constraints_ps" =>
            _require_key(model_after, "model_after", "num_linear_constraints"),
        "parity_presolve_time" =>
            _require_key(parity_stats, "parity_stats", "parity_presolve_time"),
        "residue_presolve_time" =>
            _require_key(residue_stats, "residue_stats", "residue_presolve_time"),
        "num_parity_substitutions" =>
            _require_key(parity_stats, "parity_stats", "num_parity_substitutions"),
        "num_pattern_substitutions" =>
            _require_key(parity_stats, "parity_stats", "num_pattern_substitutions"),
        "num_variables_fixed" =>
            _require_key(parity_stats, "parity_stats", "num_variables_fixed_after_parity_substitution"),
        "num_constraints_tightened" =>
            _require_key(residue_stats, "residue_stats", "num_constraints_tightened"),
        "num_constraints_eq_tightened" =>
            _require_key(residue_stats, "residue_stats", "num_constraints_eq_tightened"),
        "avg_relative_bound_tightening" =>
            _require_key(residue_stats, "residue_stats", "avg_relative_bound_tightening"),
        "infeasibility_detected" => parity_infeasible || residue_infeasible,
        "infeasibility_source" => _infeasibility_source(parity_stats, residue_stats),
    )
end

function _validate_header(instances_csv::AbstractString)
    actual_header = String.(propertynames(CSV.File(instances_csv; limit = 0)))
    missing_columns = setdiff(INSTANCE_COLUMNS, actual_header)
    isempty(missing_columns) || error(
        "Missing required column(s) in $instances_csv: $(join(missing_columns, ", "))"
    )
    return nothing
end

function merge_stats(config::CliConfig)
    isfile(config.instances_csv) || error("Instances CSV not found: $(config.instances_csv)")
    isdir(config.stats_dir) || error("Stats directory not found: $(config.stats_dir)")
    _validate_header(config.instances_csv)

    output_rows = NamedTuple[]
    seen_instance_names = Set{String}()

    for (row_index, row) in enumerate(CSV.File(config.instances_csv))
        instance_name_value = _row_value(row, "instance_name")
        instance_name_value !== missing || error(
            "Missing instance_name on row $(row_index + 1) in $(config.instances_csv)"
        )
        instance_name = string(instance_name_value)
        !isempty(instance_name) || error(
            "Empty instance_name on row $(row_index + 1) in $(config.instances_csv)"
        )
        !(instance_name in seen_instance_names) || error(
            "Duplicate instance_name in $(config.instances_csv): $instance_name"
        )
        push!(seen_instance_names, instance_name)

        stats_path = _stats_json_path(config.stats_dir, instance_name)
        isfile(stats_path) || error("Missing stats JSON for $instance_name: $stats_path")

        values = Dict{String, Any}(
            column => _row_value(row, column) for column in INSTANCE_COLUMNS
        )
        merge!(values, _stats_values(JSON.parsefile(stats_path), instance_name))
        push!(output_rows, _ordered_named_tuple(values, OUTPUT_COLUMNS))
    end

    mkpath(dirname(config.output))
    CSV.write(config.output, output_rows; writeheader = true)
    return config.output
end

function main(args::Vector{String} = copy(ARGS))
    config = parse_args(args)
    config === nothing && return nothing
    output = merge_stats(config)
    println("Wrote merged presolve stats CSV to $output")
    return output
end

end # module

if _RUNNING_AS_SCRIPT
    MergeGraphEmbeddingPresolveStats.main(copy(ARGS))
end
