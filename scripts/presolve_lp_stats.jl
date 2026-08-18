#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module PresolveLpStatsScript

using Dates
import JSON
import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC

Base.@kwdef struct CliConfig
    input_path::String
    output_path::Union{Nothing, String} = nothing
    out_dir::Union{Nothing, String} = nothing
    residue_strategy::Symbol = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY
    residue_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD
    treewidth_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD
end

function usage()
    return """
    Usage:
      julia --project=. scripts/presolve_lp_stats.jl instance.lp [options]

    Options:
      --output path                Explicit JSON output path
      --out-dir dir                Directory for the default <instance>_stats.json
      --residue-strategy name      small_primes, powers_of_two, or divisor_free
      --residue-threshold n        Maximum residue modulus threshold
      --treewidth-threshold n      Residue DP treewidth threshold
      -h, --help                   Show this help
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

function _parse_int(value::AbstractString, name::AbstractString)
    try
        parsed = parse(Int, value)
        parsed >= 0 || error()
        return parsed
    catch
        error("Invalid nonnegative Int for $name: $value")
    end
end

function _parse_strategy(value::AbstractString)
    normalized = startswith(value, ":") ? value[2:end] : value
    return Symbol(normalized)
end

function parse_args(args::Vector{String})
    if !isempty(args) && args[1] in ("-h", "--help")
        println(usage())
        return nothing
    end

    positionals = String[]
    output_path = nothing
    out_dir = nothing
    residue_strategy = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY
    residue_threshold = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD
    treewidth_threshold = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD

    index = 1
    while index <= length(args)
        arg = args[index]
        if startswith(arg, "--")
            option, value, consumed = _option_value(args, index, arg)
            if option == "--output"
                output_path = value
            elseif option == "--out-dir"
                out_dir = value
            elseif option == "--residue-strategy"
                residue_strategy = _parse_strategy(value)
            elseif option == "--residue-threshold"
                residue_threshold = _parse_int(value, option)
            elseif option == "--treewidth-threshold"
                treewidth_threshold = _parse_int(value, option)
            else
                error("Unknown option: $option")
            end
            index += consumed
        else
            push!(positionals, arg)
            index += 1
        end
    end

    length(positionals) == 1 || error("Expected exactly one LP file path, got $(length(positionals))")
    output_path === nothing || out_dir === nothing || error("--output and --out-dir are mutually exclusive")

    input_path = abspath(positionals[1])
    return CliConfig(
        input_path = input_path,
        output_path = output_path === nothing ? nothing : abspath(output_path),
        out_dir = out_dir === nothing ? nothing : abspath(out_dir),
        residue_strategy = residue_strategy,
        residue_threshold = residue_threshold,
        treewidth_threshold = treewidth_threshold,
    )
end

function default_output_path(config::CliConfig)
    config.output_path !== nothing && return config.output_path
    instance_name = splitext(basename(config.input_path))[1]
    filename = "$(instance_name)_stats.json"
    return config.out_dir === nothing ? abspath(filename) : joinpath(config.out_dir, filename)
end

function json_safe(value)
    if value === nothing || value isa Bool || value isa AbstractString
        return value
    elseif value isa AbstractFloat
        isfinite(value) && return value
        isnan(value) && return "NaN"
        return value > 0 ? "Inf" : "-Inf"
    elseif value isa Integer
        return value
    elseif value isa Symbol || value isa Enum || value isa VersionNumber
        return string(value)
    elseif value isa AbstractDict
        result = Dict{String, Any}()
        for (key, inner) in value
            result[string(key)] = json_safe(inner)
        end
        return result
    elseif value isa NamedTuple
        return Dict{String, Any}(String(key) => json_safe(getfield(value, key)) for key in keys(value))
    elseif value isa Tuple || value isa AbstractVector || value isa AbstractSet
        return Any[json_safe(inner) for inner in value]
    else
        return string(value)
    end
end

function struct_stats_dict(stats)
    return Dict{String, Any}(
        String(field) => json_safe(getfield(stats, field)) for field in fieldnames(typeof(stats))
    )
end

function _utc_timestamp()
    return string(Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss"), "Z")
end

function _build_core_model(moi_model)
    try
        return QIP.build_model(QIP.from_moi(moi_model))
    catch err
        message = sprint(showerror, err)
        if occursin("continous", message) || occursin("continuous", message)
            error(
                "Unsupported LP: QIPresolve supports only integer and binary variables; " *
                "continuous variables were found or implied. Original error: $message",
            )
        end
        rethrow()
    end
end

function _read_lp_core_model(file_path::AbstractString)
    moi_model = FF.Model(format = FF.FORMAT_LP)
    MOI.read_from_file(moi_model, file_path)
    return _build_core_model(moi_model)
end

function _legacy_ranged_constraint_match(line::AbstractString)
    return match(
        r"^(\s*)([^:\s][^:]*):\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*<=\s*(.+)\s*<=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$",
        line,
    )
end

function repair_legacy_ranged_lp_constraints(input_path::AbstractString, output_path::AbstractString)
    repaired_constraints = 0
    open(output_path, "w") do output
        for line in eachline(input_path)
            matched = _legacy_ranged_constraint_match(line)
            if matched === nothing
                println(output, line)
                continue
            end

            indent, name, lower, expr, upper = matched.captures
            clean_name = strip(name)
            clean_expr = strip(expr)
            println(output, "$(indent)$(clean_name)_lb: $(clean_expr) >= $(lower)")
            println(output, "$(indent)$(clean_name)_ub: $(clean_expr) <= $(upper)")
            repaired_constraints += 1
        end
    end
    return repaired_constraints
end

function load_lp_core_model(file_path::AbstractString)
    try
        model = _read_lp_core_model(file_path)
        return (
            model = model,
            legacy_lp_repaired = false,
            legacy_lp_repaired_constraints = 0,
            legacy_lp_repair_error = nothing,
        )
    catch initial_error
        repair_path, io = mktemp()
        close(io)
        try
            repaired_constraints = repair_legacy_ranged_lp_constraints(file_path, repair_path)
            repaired_constraints > 0 || throw(initial_error)

            try
                model = _read_lp_core_model(repair_path)
                return (
                    model = model,
                    legacy_lp_repaired = true,
                    legacy_lp_repaired_constraints = repaired_constraints,
                    legacy_lp_repair_error = sprint(showerror, initial_error),
                )
            catch repair_error
                error(
                    "Failed to read LP file. Initial error: $(sprint(showerror, initial_error)). " *
                    "Legacy ranged-constraint repair also failed: $(sprint(showerror, repair_error))",
                )
            end
        catch repair_error
            repair_error === initial_error && throw(initial_error)
            rethrow()
        finally
            rm(repair_path; force = true)
        end
    end
end

function run_tracked_presolve!(model::PC.QPModel, config::CliConfig)
    postsolver = PC.ParityPostsolver(keys(model.vars))
    parity_stats_accumulator = PC._ParityStatsAccumulator()
    residue_stats_accumulator = PC._ResidueStatsAccumulator()

    PC.normalize!(model, postsolver)
    model_before = PC.ModelStats(model)

    if !model.infeasible
        moduli = PC._generate_residue_moduli(config.residue_strategy, config.residue_threshold)
        propagator = PC.PropagationManager(PC.VarId[])

        PC.parity_presolve!(model, propagator, postsolver, parity_stats_accumulator)

        if !model.infeasible
            residue_result = PC._residue_presolve_pass!(
                model,
                moduli,
                nothing,
                postsolver,
                residue_stats_accumulator;
                treewidth_threshold = config.treewidth_threshold,
            )

            while !model.infeasible && residue_result.tightened_to_equality
                parity_result = PC.parity_presolve!(
                    model,
                    propagator,
                    postsolver,
                    parity_stats_accumulator,
                )
                model.infeasible && break
                parity_result.domains_changed || break

                residue_result = PC._residue_presolve_pass!(
                    model,
                    moduli,
                    parity_result.coefficient_changed_constraint_ids,
                    postsolver,
                    residue_stats_accumulator;
                    treewidth_threshold = config.treewidth_threshold,
                )
            end
        end
    end

    return (
        model = model,
        postsolver = postsolver,
        model_before = model_before,
        model_after = PC.ModelStats(model),
        parity_stats = parity_stats_accumulator.stats,
        residue_stats = residue_stats_accumulator.stats,
    )
end

function parameters_dict(config::CliConfig)
    return Dict{String, Any}(
        "residue_strategy" => string(config.residue_strategy),
        "residue_threshold" => config.residue_threshold,
        "treewidth_threshold" => config.treewidth_threshold,
    )
end

function metadata_dict(config::CliConfig, output_path::AbstractString, load_info)
    return Dict{String, Any}(
        "generated_at_utc" => _utc_timestamp(),
        "input_path" => config.input_path,
        "input_size_bytes" => filesize(config.input_path),
        "output_path" => output_path,
        "julia_version" => string(VERSION),
        "script" => basename(@__FILE__),
        "legacy_lp_repaired" => load_info.legacy_lp_repaired,
        "legacy_lp_repaired_constraints" => load_info.legacy_lp_repaired_constraints,
        "legacy_lp_repair_error" => load_info.legacy_lp_repair_error,
    )
end

function write_json(path::AbstractString, data)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, json_safe(data), 4)
        println(io)
    end
    return path
end

function run_analysis(config::CliConfig)
    isfile(config.input_path) || error("LP file not found: $(config.input_path)")

    output_path = default_output_path(config)
    instance_name = splitext(basename(config.input_path))[1]

    load_info = load_lp_core_model(config.input_path)
    loaded_model = load_info.model
    model_loaded = PC.ModelStats(loaded_model)

    presolve_model = deepcopy(loaded_model)
    presolve_result = run_tracked_presolve!(presolve_model, config)

    data = Dict{String, Any}(
        "instance_name" => instance_name,
        "metadata" => metadata_dict(config, output_path, load_info),
        "parameters" => parameters_dict(config),
        "model_loaded" => struct_stats_dict(model_loaded),
        "model_before" => struct_stats_dict(presolve_result.model_before),
        "model_after" => struct_stats_dict(presolve_result.model_after),
        "parity_stats" => struct_stats_dict(presolve_result.parity_stats),
        "residue_stats" => struct_stats_dict(presolve_result.residue_stats),
    )

    write_json(output_path, data)
    return (output_path = output_path, data = data)
end

function main(args::Vector{String} = copy(ARGS))
    config = parse_args(args)
    config === nothing && return nothing
    result = run_analysis(config)
    println("Wrote stats JSON to $(result.output_path)")
    return result.output_path
end

end # module

if _RUNNING_AS_SCRIPT
    PresolveLpStatsScript.main(copy(ARGS))
end
