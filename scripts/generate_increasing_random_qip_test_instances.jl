#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using Dates
import JuMP
import MathOptInterface as MOI

using QIPresolve.InstanceGeneration: generate_random_qip_model
using QIPresolve.ModelIO: save_moi

const DEFAULT_TARGET = joinpath("instances", "test")
const DEFAULT_SEED_BASE = 10_000
const DEFAULT_MAX_ATTEMPTS = 100
const INSTANCE_PREFIX = "test_qip"
const NVAR_VALUES = collect(10:10:100)
const CONSTRAINT_SLACK_RANGE_LABEL = "-10:10"

const CSV_COLUMNS = [
    "instance_name",
    "file_name",
    "num",
    "created_at",
    "nvars",
    "ncons",
    "seed",
    "eq_constraints",
    "ineq_constraints",
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
    "target" => :target,
    "csv" => :csv,
    "seed-base" => :seed_base,
    "max-attempts" => :max_attempts,
)

const RANDOM_QIP_KWARGS = (
    p_con_eq = 0.5,
    var_threshold_lb = -10,
    var_threshold_ub = 10,
    p_var_is_candidate = 0.02,
    p_var_bilin = 0.4,
    p_var_diag = 0.5,
    p_var_lin = 0.0,
    coeff_lb = -50,
    coeff_ub = 50,
    force_diag_even = false,
    force_lin_even = false,
    force_feasibility = true,
    constraint_slack_range = collect(-10:10),
)

Base.@kwdef struct GeneratorConfig
    target::String = abspath(DEFAULT_TARGET)
    csv_path::String = joinpath(abspath(DEFAULT_TARGET), "instances.csv")
    seed_base::Int = DEFAULT_SEED_BASE
    max_attempts::Int = DEFAULT_MAX_ATTEMPTS
end

struct GeneratedInstance
    instance_name::String
    file_name::String
    file_path::String
    num::Int
    nvars::Int
    ncons::Int
    seed::Int
    eq_constraints::Int
    ineq_constraints::Int
    model::JuMP.Model
end

function usage()
    println("Usage:")
    println("  julia --project=. scripts/generate_increasing_random_qip_test_instances.jl [options]")
    println()
    println("Options:")
    println("  --target $(DEFAULT_TARGET)")
    println("  --csv <target>/instances.csv")
    println("  --seed-base $(DEFAULT_SEED_BASE)")
    return println("  --max-attempts $(DEFAULT_MAX_ATTEMPTS)")
end

function parse_int(value::AbstractString, name::AbstractString)::Int
    try
        return parse(Int, value)
    catch
        error("Invalid $name: $value")
    end
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

    config = GeneratorConfig(
        target = target,
        csv_path = csv_path,
        seed_base = parse_int(get(options, :seed_base, string(DEFAULT_SEED_BASE)), "seed_base"),
        max_attempts = parse_int(get(options, :max_attempts, string(DEFAULT_MAX_ATTEMPTS)), "max_attempts"),
    )

    config.max_attempts >= 1 || error("max_attempts must be >= 1")
    return config
end

instance_name(nvars::Int) = "$(INSTANCE_PREFIX)_n$nvars"

function planned_file_paths(config::GeneratorConfig)
    return [
        joinpath(config.target, "$(instance_name(nvars)).lp")
        for nvars in NVAR_VALUES
    ]
end

function ensure_output_is_new!(config::GeneratorConfig)
    if isfile(config.csv_path) && filesize(config.csv_path) > 0
        error("Refusing to overwrite existing nonempty CSV: $(config.csv_path)")
    end

    for file_path in planned_file_paths(config)
        isfile(file_path) && error("Refusing to overwrite existing instance file: $file_path")
    end
    return nothing
end

function constraint_counts(model::JuMP.Model)
    eq_count = length(JuMP.all_constraints(model, JuMP.QuadExpr, MOI.EqualTo{Float64}))
    ineq_count = length(JuMP.all_constraints(model, JuMP.QuadExpr, MOI.Interval{Float64}))
    return eq_count, ineq_count
end

function seed_for_attempt(config::GeneratorConfig, instance_idx::Int, attempt_idx::Int)
    return config.seed_base + (instance_idx - 1) * config.max_attempts + attempt_idx
end

function generate_instance(config::GeneratorConfig, instance_idx::Int, nvars::Int)
    ncons = 2 * nvars
    for attempt_idx in 0:(config.max_attempts - 1)
        seed = seed_for_attempt(config, instance_idx, attempt_idx)
        model, _ = generate_random_qip_model(
            nvars,
            ncons;
            RANDOM_QIP_KWARGS...,
            seed = seed,
        )

        eq_count, ineq_count = constraint_counts(model)
        if eq_count > 0 && ineq_count > 0
            name = instance_name(nvars)
            file_name = "$name.lp"
            return GeneratedInstance(
                name,
                file_name,
                joinpath(config.target, file_name),
                instance_idx,
                nvars,
                ncons,
                seed,
                eq_count,
                ineq_count,
                model,
            )
        end
    end

    error("Failed to generate nvars=$nvars with both equality and inequality constraints after $(config.max_attempts) attempts")
end

function blank_csv_row()
    return Dict(column => "" for column in CSV_COLUMNS)
end

function csv_row(instance::GeneratedInstance)
    row = blank_csv_row()
    row["instance_name"] = instance.instance_name
    row["file_name"] = instance.file_name
    row["num"] = string(instance.num)
    row["created_at"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    row["nvars"] = string(instance.nvars)
    row["ncons"] = string(instance.ncons)
    row["seed"] = string(instance.seed)
    row["eq_constraints"] = string(instance.eq_constraints)
    row["ineq_constraints"] = string(instance.ineq_constraints)
    row["p_con_eq"] = string(RANDOM_QIP_KWARGS.p_con_eq)
    row["var_threshold_lb"] = string(RANDOM_QIP_KWARGS.var_threshold_lb)
    row["var_threshold_ub"] = string(RANDOM_QIP_KWARGS.var_threshold_ub)
    row["p_var_is_candidate"] = string(RANDOM_QIP_KWARGS.p_var_is_candidate)
    row["p_var_bilin"] = string(RANDOM_QIP_KWARGS.p_var_bilin)
    row["p_var_diag"] = string(RANDOM_QIP_KWARGS.p_var_diag)
    row["p_var_lin"] = string(RANDOM_QIP_KWARGS.p_var_lin)
    row["coeff_lb"] = string(RANDOM_QIP_KWARGS.coeff_lb)
    row["coeff_ub"] = string(RANDOM_QIP_KWARGS.coeff_ub)
    row["force_diag_even"] = string(RANDOM_QIP_KWARGS.force_diag_even)
    row["force_lin_even"] = string(RANDOM_QIP_KWARGS.force_lin_even)
    row["force_feasibility"] = string(RANDOM_QIP_KWARGS.force_feasibility)
    row["constraint_slack_range"] = CONSTRAINT_SLACK_RANGE_LABEL
    return (; (Symbol(column) => row[column] for column in CSV_COLUMNS)...)
end

function run(config::GeneratorConfig)
    mkpath(config.target)
    mkpath(dirname(config.csv_path))
    ensure_output_is_new!(config)

    generated_instances = [
        generate_instance(config, idx, nvars)
        for (idx, nvars) in enumerate(NVAR_VALUES)
    ]

    for instance in generated_instances
        save_moi(JuMP.backend(instance.model), instance.file_path)
        println("saved $(instance.instance_name) -> $(instance.file_path)")
    end

    CSV.write(config.csv_path, csv_row.(generated_instances); writeheader = true)
    println("updated csv: $(config.csv_path)")
    return nothing
end

function main(args::Vector{String})
    config = build_config(args)
    return run(config)
end

main(ARGS)
