#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using Dates
using JuMP: backend

using QIPresolve.InstanceGeneration: generate_random_qip_model
using QIPresolve.ModelIO: save_moi

const DEFAULT_TARGET = joinpath("results", "random_qip_instances")
const INSTANCE_PREFIX = "qip"

const CSV_COLUMNS = [
    "instance_name",
    "num",
    "created_at",
    "nvars",
    "ncons",
    "seed",
    "p_con_eq",
    "var_threshold_lb",
    "var_threshold_ub",
    "p_var_is_candidate",
    "p_var_bilin",
    "p_var_diag",
    "p_var_lin",
    "coeff_lb",
    "coeff_ub",
    "force_diag_even",
    "force_lin_even",
    "force_feasibility",
    "constraint_slack_range",
]

const CLI_KEYS = Dict(
    "count" => :count,
    "target" => :target,
    "csv" => :csv,
    "instance-prefix" => :instance_name_prefix,
    "instance-name-prefix" => :instance_name_prefix,
    "name-prefix" => :instance_name_prefix,
    "nvars" => :nvars,
    "n-vars" => :nvars,
    "ncons" => :ncons,
    "n-cons" => :ncons,
    "nconstraints" => :ncons,
    "n-constraints" => :ncons,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "p-con-eq" => :p_con_eq,
    "var-threshold-lb" => :var_threshold_lb,
    "var-threshold-ub" => :var_threshold_ub,
    "p-var-is-candidate" => :p_var_is_candidate,
    "p-var-candidate" => :p_var_is_candidate,
    "p-var-candiate" => :p_var_is_candidate,
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

Base.@kwdef struct GeneratorConfig
    count::Int = 1
    target::String = abspath(DEFAULT_TARGET)
    csv_path::String = joinpath(abspath(DEFAULT_TARGET), "instances.csv")
    instance_name_prefix::String = ""
    nvars::Int = 100
    ncons::Int = 200
    seed_base::Int = 1
    seed_step::Int = 1
    p_con_eq::Float64 = 0.0
    var_threshold_lb::Int = -10
    var_threshold_ub::Int = 10
    p_var_is_candidate::Float64 = 0.02
    p_var_bilin::Float64 = 0.4
    p_var_diag::Float64 = 0.5
    p_var_lin::Float64 = 0.0
    coeff_lb::Int = -50
    coeff_ub::Int = 50
    force_diag_even::Bool = false
    force_lin_even::Bool = false
    force_feasibility::Bool = true
    constraint_slack_range::Vector{Int} = collect(-10:10)
    constraint_slack_range_label::String = "-10:10"
end

function usage()
    println("Usage:")
    println("  julia --project=. scripts/generate_random_qip_instances.jl [options]")
    println()
    println("Common options:")
    println("  --count 1")
    println("  --target $(DEFAULT_TARGET)")
    println("  --csv <target>/instances.csv")
    println("  --instance-prefix <prefix>")
    println("  --nvars 100")
    println("  --ncons 200")
    println("  --seed-base 1")
    println("  --seed-step 1")
    println()
    println("Random QIP options:")
    println("  --p-con-eq 0.0")
    println("  --var-threshold-lb -10")
    println("  --var-threshold-ub 10")
    println("  --p-var-is-candidate 0.02")
    println("  --p-var-bilin 0.4")
    println("  --p-var-diag 0.5")
    println("  --p-var-lin 0.0")
    println("  --coeff-lb -50")
    println("  --coeff-ub 50")
    println("  --force-diag-even false")
    println("  --force-lin-even false")
    println("  --force-feasibility true")
    return println("  --constraint-slack-range -10:10")
end

function parse_int(value::AbstractString, name::AbstractString)::Int
    try
        return parse(Int, value)
    catch
        error("Invalid $name: $value")
    end
end

function parse_float(value::AbstractString, name::AbstractString)::Float64
    try
        return parse(Float64, value)
    catch
        error("Invalid $name: $value")
    end
end

function parse_bool(value::AbstractString, name::AbstractString)::Bool
    normalized = lowercase(strip(value))
    normalized in ("true", "1", "yes") && return true
    normalized in ("false", "0", "no") && return false
    error("Invalid $name: $value. Expected true/false, 1/0, or yes/no.")
end

function parse_instance_name_prefix(value::AbstractString)::String
    prefix = strip(value)
    if occursin("/", prefix) || occursin("\\", prefix)
        error("Invalid --instance-prefix: path separators are not allowed")
    end
    return prefix
end

function parse_constraint_slack_range(value::AbstractString)
    normalized = replace(strip(value), " " => "")
    isempty(normalized) && error("Invalid constraint_slack_range: value must be nonempty")

    if occursin(":", normalized)
        parts = split(normalized, ":"; limit = 2)
        length(parts) == 2 || error("Invalid constraint_slack_range: $value")
        lower = parse_int(parts[1], "constraint_slack_range lower bound")
        upper = parse_int(parts[2], "constraint_slack_range upper bound")
        lower <= upper || error("constraint_slack_range lower bound must be <= upper bound")
        return collect(lower:upper), "$lower:$upper"
    end

    values = [parse_int(part, "constraint_slack_range value") for part in split(normalized, ",")]
    isempty(values) && error("Invalid constraint_slack_range: value must be nonempty")
    return values, join(values, ",")
end

function parse_raw_options(args::Vector{String})
    options = Dict{Symbol, String}()
    idx = 1
    while idx <= length(args)
        arg = args[idx]
        if arg in ("-h", "--help")
            usage()
            exit(0)
        end
        startswith(arg, "--") || error("Unexpected positional argument: $arg")

        raw_key_value = arg[3:end]
        raw_key, value = if occursin("=", raw_key_value)
            split(raw_key_value, "="; limit = 2)
        else
            idx += 1
            idx <= length(args) || error("Missing value for $arg")
            raw_key_value, args[idx]
        end

        lookup_key = lowercase(replace(raw_key, "_" => "-"))
        key = get(CLI_KEYS, lookup_key, nothing)
        key === nothing && error("Unknown option: --$raw_key")
        options[key] = value
        idx += 1
    end
    return options
end

function build_config(args::Vector{String})::GeneratorConfig
    options = parse_raw_options(args)
    target = abspath(get(options, :target, DEFAULT_TARGET))
    csv_path = haskey(options, :csv) ?
        abspath(options[:csv]) :
        joinpath(target, "instances.csv")
    constraint_slack_range, constraint_slack_range_label = parse_constraint_slack_range(
        get(options, :constraint_slack_range, "-10:10")
    )

    config = GeneratorConfig(
        count = parse_int(get(options, :count, "1"), "count"),
        target = target,
        csv_path = csv_path,
        instance_name_prefix = parse_instance_name_prefix(get(options, :instance_name_prefix, "")),
        nvars = parse_int(get(options, :nvars, "100"), "nvars"),
        ncons = parse_int(get(options, :ncons, "200"), "ncons"),
        seed_base = parse_int(get(options, :seed_base, "1"), "seed_base"),
        seed_step = parse_int(get(options, :seed_step, "1"), "seed_step"),
        p_con_eq = parse_float(get(options, :p_con_eq, "0.0"), "p_con_eq"),
        var_threshold_lb = parse_int(get(options, :var_threshold_lb, "-10"), "var_threshold_lb"),
        var_threshold_ub = parse_int(get(options, :var_threshold_ub, "10"), "var_threshold_ub"),
        p_var_is_candidate = parse_float(get(options, :p_var_is_candidate, "0.02"), "p_var_is_candidate"),
        p_var_bilin = parse_float(get(options, :p_var_bilin, "0.4"), "p_var_bilin"),
        p_var_diag = parse_float(get(options, :p_var_diag, "0.5"), "p_var_diag"),
        p_var_lin = parse_float(get(options, :p_var_lin, "0.0"), "p_var_lin"),
        coeff_lb = parse_int(get(options, :coeff_lb, "-50"), "coeff_lb"),
        coeff_ub = parse_int(get(options, :coeff_ub, "50"), "coeff_ub"),
        force_diag_even = parse_bool(get(options, :force_diag_even, "false"), "force_diag_even"),
        force_lin_even = parse_bool(get(options, :force_lin_even, "false"), "force_lin_even"),
        force_feasibility = parse_bool(get(options, :force_feasibility, "true"), "force_feasibility"),
        constraint_slack_range = constraint_slack_range,
        constraint_slack_range_label = constraint_slack_range_label,
    )

    config.count >= 1 || error("count must be >= 1")
    return config
end

function ensure_csv_header!(path::AbstractString)
    (!isfile(path) || filesize(path) == 0) && return nothing

    actual_header = String.(propertynames(CSV.File(path; limit = 0)))
    actual_header == CSV_COLUMNS || error(
        "CSV header does not match expected schema in $path"
    )
    return nothing
end

function parse_instance_num(instance_name::AbstractString, instance_name_prefix::AbstractString)
    stem = "$(instance_name_prefix)$(INSTANCE_PREFIX)_"
    startswith(instance_name, stem) || return nothing
    suffix = SubString(instance_name, ncodeunits(stem) + 1)
    return tryparse(Int, suffix)
end

function next_instance_number(csv_path::AbstractString, instance_name_prefix::AbstractString)::Int
    !isfile(csv_path) && return 1
    filesize(csv_path) == 0 && return 1

    max_num = 0
    for (row_idx, row) in enumerate(CSV.File(csv_path))
        instance_name = row.instance_name === missing ? "" : string(row.instance_name)
        name_num = parse_instance_num(instance_name, instance_name_prefix)
        name_num === nothing && continue

        parsed_num = row.num === missing ? nothing : tryparse(Int, string(row.num))
        if parsed_num === nothing
            parsed_num = name_num
        end
        parsed_num === nothing && error(
            "Invalid num for prefix $(instance_name_prefix)$INSTANCE_PREFIX on CSV row $(row_idx + 1) in $csv_path"
        )
        max_num = max(max_num, parsed_num)
    end

    return max_num + 1
end

function random_qip_kwargs(config::GeneratorConfig)
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

function blank_csv_row()
    return Dict(column => "" for column in CSV_COLUMNS)
end

function csv_row(
        config::GeneratorConfig,
        instance_name::String,
        num::Int,
        seed::Int,
    )
    row = blank_csv_row()
    row["instance_name"] = instance_name
    row["num"] = string(num)
    row["created_at"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    row["nvars"] = string(config.nvars)
    row["ncons"] = string(config.ncons)
    row["seed"] = string(seed)
    row["p_con_eq"] = string(config.p_con_eq)
    row["var_threshold_lb"] = string(config.var_threshold_lb)
    row["var_threshold_ub"] = string(config.var_threshold_ub)
    row["p_var_is_candidate"] = string(config.p_var_is_candidate)
    row["p_var_bilin"] = string(config.p_var_bilin)
    row["p_var_diag"] = string(config.p_var_diag)
    row["p_var_lin"] = string(config.p_var_lin)
    row["coeff_lb"] = string(config.coeff_lb)
    row["coeff_ub"] = string(config.coeff_ub)
    row["force_diag_even"] = string(config.force_diag_even)
    row["force_lin_even"] = string(config.force_lin_even)
    row["force_feasibility"] = string(config.force_feasibility)
    row["constraint_slack_range"] = config.constraint_slack_range_label
    return (; (Symbol(column) => row[column] for column in CSV_COLUMNS)...)
end

function append_csv_row!(path::AbstractString, row::NamedTuple)
    has_existing_rows = isfile(path) && filesize(path) > 0
    CSV.write(path, [row]; append = has_existing_rows, writeheader = !has_existing_rows)
    return nothing
end

function run(config::GeneratorConfig)
    mkpath(config.target)
    mkpath(dirname(config.csv_path))
    ensure_csv_header!(config.csv_path)

    start_num = next_instance_number(config.csv_path, config.instance_name_prefix)
    for offset in 0:(config.count - 1)
        num = start_num + offset
        seed = config.seed_base + (num - 1) * config.seed_step
        instance_name = "$(config.instance_name_prefix)$(INSTANCE_PREFIX)_$num"
        file_name = "$instance_name.lp"
        file_path = joinpath(config.target, file_name)
        isfile(file_path) && error("Refusing to overwrite existing instance file: $file_path")

        model, _ = generate_random_qip_model(
            config.nvars,
            config.ncons;
            random_qip_kwargs(config)...,
            seed = seed,
        )
        save_moi(backend(model), file_path)
        append_csv_row!(
            config.csv_path,
            csv_row(config, instance_name, num, seed),
        )
        println("saved $instance_name -> $file_path")
    end

    println("updated csv: $(config.csv_path)")
    return nothing
end

function main(args::Vector{String})
    config = build_config(args)
    return run(config)
end

main(ARGS)
