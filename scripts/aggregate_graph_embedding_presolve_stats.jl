#!/usr/bin/env julia

const _RUNNING_AGGREGATE_PRESOLVE_STATS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AGGREGATE_PRESOLVE_STATS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module AggregateGraphEmbeddingPresolveStats

import CSV
import Statistics

const DEFAULT_BASE_DIR = joinpath("results", "graph_embedding_instances")
const DEFAULT_INPUT = joinpath(DEFAULT_BASE_DIR, "instances_with_presolve_stats.csv")
const DEFAULT_OUTPUT =
    joinpath(DEFAULT_BASE_DIR, "instances_with_presolve_stats_aggregated.csv")
const OUTPUT_FLOAT_DIGITS = 3

const GROUP_COLUMNS = [
    "type",
    "n",
    "R",
    "num_anchors",
    "alpha",
    "edge_density",
    "infeas_strategy",
    "infeas_base",
    "box_scale",
]

const AVERAGE_COLUMNS = [
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
]

const INFEASIBILITY_SOURCES = ["none", "elimination", "propagation", "residue"]
const COUNTED_INFEASIBILITY_SOURCES = ["elimination", "propagation", "residue"]
const RELATIVE_LOG_DOMAIN_REDUCTION_COLUMN = "relative_log_domain_reduction_percent"
const AVG_RELATIVE_BOUND_TIGHTENING_COLUMN = "avg_relative_bound_tightening"
const AVG_RELATIVE_BOUND_TIGHTENING_PERCENT_COLUMN =
    "avg_relative_bound_tightening_percent"
const SOURCE_COUNT_COLUMNS = [
    "infeasibility_source_$(source)_count" for source in COUNTED_INFEASIBILITY_SOURCES
]
const REQUIRED_COLUMNS = [
    ["instance_name"];
    GROUP_COLUMNS;
    AVERAGE_COLUMNS;
    ["infeasibility_detected", "infeasibility_source"];
]
const OUTPUT_COLUMNS = [
    GROUP_COLUMNS;
    ["num_instances"];
    AVERAGE_COLUMNS[1:2];
    [RELATIVE_LOG_DOMAIN_REDUCTION_COLUMN];
    AVERAGE_COLUMNS[3:(end - 1)];
    [AVG_RELATIVE_BOUND_TIGHTENING_PERCENT_COLUMN];
    ["infeasibility_detected_count"];
    SOURCE_COUNT_COLUMNS;
]

Base.@kwdef struct CliConfig
    input::String = abspath(DEFAULT_INPUT)
    output::String = abspath(DEFAULT_OUTPUT)
end

function usage()
    return """
    Usage:
      julia --project=. scripts/aggregate_graph_embedding_presolve_stats.jl [options]

    Options:
      --input path     Input merged graph embedding presolve stats CSV
      --output path    Output aggregate CSV path
      -h, --help       Show this help
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
    input = DEFAULT_INPUT
    output = DEFAULT_OUTPUT

    index = 1
    while index <= length(args)
        arg = args[index]
        if arg in ("-h", "--help")
            println(usage())
            return nothing
        end
        startswith(arg, "--") || error("Unexpected positional argument: $arg")

        option, value, consumed = _option_value(args, index, arg)
        if option == "--input"
            input = value
        elseif option == "--output"
            output = value
        else
            error("Unknown option: $option")
        end
        index += consumed
    end

    return CliConfig(input = abspath(input), output = abspath(output))
end

function _row_value(row, column::AbstractString)
    return getproperty(row, Symbol(column))
end

function _ordered_named_tuple(values::AbstractDict, columns::Vector{String})
    names = Tuple(Symbol.(columns))
    rounded_values = Tuple(
        value isa AbstractFloat ? round(value; digits = OUTPUT_FLOAT_DIGITS) : value
        for value in (values[column] for column in columns)
    )
    return NamedTuple{names}(rounded_values)
end

function _validate_header(input::AbstractString)
    actual_header = String.(propertynames(CSV.File(input; limit = 0)))
    missing_columns = setdiff(REQUIRED_COLUMNS, actual_header)
    isempty(missing_columns) || error(
        "Missing required column(s) in $input: $(join(missing_columns, ", "))"
    )
    return nothing
end

function _row_context(row, row_index::Int, input::AbstractString)
    instance_name = _row_value(row, "instance_name")
    instance_name === missing && error(
        "Missing instance_name on row $(row_index + 1) in $input"
    )
    name = string(instance_name)
    isempty(name) && error("Empty instance_name on row $(row_index + 1) in $input")
    return "instance $name on row $(row_index + 1) in $input"
end

function _validate_row(row, row_index::Int, input::AbstractString)
    context = _row_context(row, row_index, input)

    for column in AVERAGE_COLUMNS
        value = _row_value(row, column)
        value isa Real && !(value isa Bool) || error(
            "Expected numeric value for $column for $context, found $(repr(value))"
        )
        isfinite(Float64(value)) || error(
            "Expected finite numeric value for $column for $context, found $(repr(value))"
        )
    end

    detected = _row_value(row, "infeasibility_detected")
    detected isa Bool || error(
        "Expected Boolean infeasibility_detected for $context, found $(repr(detected))"
    )

    source_value = _row_value(row, "infeasibility_source")
    source_value isa AbstractString || error(
        "Expected string infeasibility_source for $context, found $(repr(source_value))"
    )
    source = string(source_value)
    source in INFEASIBILITY_SOURCES || error(
        "Unknown infeasibility_source for $context: $(repr(source))"
    )
    detected == (source != "none") || error(
        "Inconsistent infeasibility fields for $context: " *
        "infeasibility_detected=$detected, infeasibility_source=$(repr(source))"
    )
    return nothing
end

function _group_key(row)
    return Tuple(_row_value(row, column) for column in GROUP_COLUMNS)
end

function aggregate_stats(config::CliConfig)
    isfile(config.input) || error("Input CSV not found: $(config.input)")
    _validate_header(config.input)

    rows = collect(CSV.File(config.input))
    isempty(rows) && error("Input CSV contains no data rows: $(config.input)")

    group_indices = Dict{Tuple, Vector{Int}}()
    ordered_keys = Tuple[]
    for (row_index, row) in enumerate(rows)
        _validate_row(row, row_index, config.input)
        key = _group_key(row)
        if !haskey(group_indices, key)
            group_indices[key] = Int[]
            push!(ordered_keys, key)
        end
        push!(group_indices[key], row_index)
    end

    output_rows = NamedTuple[]
    for key in ordered_keys
        indices = group_indices[key]
        first_row = rows[first(indices)]
        values = Dict{String, Any}(
            column => _row_value(first_row, column) for column in GROUP_COLUMNS
        )
        values["num_instances"] = length(indices)

        for column in AVERAGE_COLUMNS
            values[column] = Statistics.mean(
                Float64(_row_value(rows[index], column)) for index in indices
            )
        end
        original_log_domain = values["log_domain_sum_orig"]
        iszero(original_log_domain) && error(
            "Cannot compute $RELATIVE_LOG_DOMAIN_REDUCTION_COLUMN for group $(repr(key)): " *
            "average log_domain_sum_orig is zero"
        )
        values[RELATIVE_LOG_DOMAIN_REDUCTION_COLUMN] =
            100.0 * (original_log_domain - values["log_domain_sum_ps"]) / original_log_domain
        values[AVG_RELATIVE_BOUND_TIGHTENING_PERCENT_COLUMN] =
            100.0 * values[AVG_RELATIVE_BOUND_TIGHTENING_COLUMN]

        values["infeasibility_detected_count"] = count(
            index -> _row_value(rows[index], "infeasibility_detected"),
            indices,
        )
        for source in COUNTED_INFEASIBILITY_SOURCES
            values["infeasibility_source_$(source)_count"] = count(
                index -> string(_row_value(rows[index], "infeasibility_source")) == source,
                indices,
            )
        end
        push!(output_rows, _ordered_named_tuple(values, OUTPUT_COLUMNS))
    end

    mkpath(dirname(config.output))
    CSV.write(config.output, output_rows; writeheader = true)
    return config.output
end

function main(args::Vector{String} = copy(ARGS))
    config = parse_args(args)
    config === nothing && return nothing
    output = aggregate_stats(config)
    println("Wrote aggregate presolve stats CSV to $output")
    return output
end

end # module

if _RUNNING_AGGREGATE_PRESOLVE_STATS_SCRIPT
    AggregateGraphEmbeddingPresolveStats.main(copy(ARGS))
end
