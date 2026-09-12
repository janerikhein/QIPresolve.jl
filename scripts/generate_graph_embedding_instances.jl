#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using Dates
using JuMP: backend

using QIPresolve.InstanceGeneration:
    BoundingBoxCandidateRejected,
    generate_2_connected_instance,
    generate_globally_rigid_instance,
    generate_laman_instance,
    generate_likely_infeasible_embedding_instance
using QIPresolve.ModelIO: save_moi

const DEFAULT_TARGET = joinpath("results", "graph_embedding_instances")
const MAX_BOUNDING_BOX_CANDIDATE_TRIES = 10_000

const CSV_COLUMNS = [
    "instance_name",
    "type",
    "num",
    "created_at",
    "n",
    "R",
    "seed",
    "num_anchors",
    "alpha",
    "edge_density",
    "pH2",
    "max_coord_tries",
    "max_tries_H2",
    "infeas_strategy",
    "infeas_base",
    "box_margin",
]

const CLI_KEYS = Dict(
    "type" => :type,
    "count" => :count,
    "target" => :target,
    "csv" => :csv,
    "instance-prefix" => :instance_name_prefix,
    "instance-name-prefix" => :instance_name_prefix,
    "name-prefix" => :instance_name_prefix,
    "n" => :n,
    "r" => :R,
    "num-anchors" => :num_anchors,
    "alpha" => :alpha,
    "seed-base" => :seed_base,
    "seed-step" => :seed_step,
    "edge-density" => :edge_density,
    "ph2" => :pH2,
    "max-coord-tries" => :max_coord_tries,
    "max-global-tries" => :max_global_tries,
    "max-tries-h1" => :max_tries_H1,
    "max-tries-h2" => :max_tries_H2,
    "infeas-strategy" => :infeas_strategy,
    "infeas-base" => :infeas_base,
    "box-margin" => :box_margin,
    "contraction-vertices" => :contraction_vertices,
)

Base.@kwdef struct GeneratorConfig
    prefix::String
    type_label::String
    count::Int = 1
    target::String = abspath(DEFAULT_TARGET)
    csv_path::String = joinpath(abspath(DEFAULT_TARGET), "instances.csv")
    instance_name_prefix::String = ""
    n::Int = 10
    R::Int = 100
    num_anchors::Int = 0
    alpha::Float64 = 0.0
    seed_base::Int = 1
    seed_step::Int = 1
    edge_density::Float64 = 0.5
    pH2::Float64 = 0.5
    max_coord_tries::Int = 10_000
    max_global_tries::Int = 10_000
    max_tries_H1::Int = 200
    max_tries_H2::Int = 300
    infeas_strategy::Symbol = :bounding_box
    infeas_base::Symbol = :globally_rigid
    box_margin::Int = 1
    contraction_vertices::Union{Nothing, NTuple{2, Int}} = nothing
    contraction_vertices_label::String = "auto"
end

function usage()
    println("Usage:")
    println("  julia --project=. scripts/generate_graph_embedding_instances.jl --type con|lam|gr|infeas [options]")
    println()
    println("Required:")
    println("  --type con|lam|gr|infeas")
    println()
    println("Common options:")
    println("  --count 1")
    println("  --target $(DEFAULT_TARGET)")
    println("  --csv <target>/instances.csv")
    println("  --instance-prefix <prefix>")
    println("  --n 10")
    println("  --R 100")
    println("  --num-anchors 0")
    println("  --alpha 0.0")
    println("  --seed-base 1")
    println("  --seed-step 1")
    println()
    println("Type-specific options:")
    println("  --edge-density 0.5")
    println("  --pH2 0.5")
    println("  --max-coord-tries 10000")
    println("  --max-global-tries 10000")
    println("  --max-tries-H1 200")
    println("  --max-tries-H2 300")
    println("  --infeas-strategy bounding_box|vertex_contraction")
    println("  --infeas-base globally_rigid|laman|random_2_connected")
    println("  --box-margin 1")
    return println("  --contraction-vertices auto|u,v")
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

function parse_instance_name_prefix(value::AbstractString)::String
    prefix = strip(value)
    if occursin("/", prefix) || occursin("\\", prefix)
        error("Invalid --instance-prefix: path separators are not allowed")
    end
    return prefix
end

function normalize_type(value::AbstractString)::Tuple{String, String}
    normalized = lowercase(replace(strip(value), "_" => "-"))
    if normalized in ("con", "2-connected", "two-connected")
        return "con", "2_connected"
    elseif normalized in ("lam", "laman")
        return "lam", "laman"
    elseif normalized in ("gr", "globally-rigid", "global-rigid")
        return "gr", "globally_rigid"
    elseif normalized in ("infeas", "infeasible")
        return "infeas", "infeasible"
    end
    error("Invalid --type: $value. Expected con, lam, gr, or infeas.")
end

function normalize_infeas_strategy(value::AbstractString)::Symbol
    normalized = lowercase(replace(strip(value), "-" => "_"))
    if normalized == "bounding_box"
        return :bounding_box
    elseif normalized in ("vertex_contraction", "contraction")
        return :vertex_contraction
    end
    error("Invalid --infeas-strategy: $value. Expected bounding_box or vertex_contraction.")
end

function normalize_infeas_base(value::AbstractString)::Symbol
    normalized = lowercase(replace(strip(value), "-" => "_"))
    if normalized in ("gr", "globally_rigid", "global_rigid")
        return :globally_rigid
    elseif normalized in ("lam", "laman")
        return :laman
    elseif normalized in ("con", "random_2_connected", "two_connected", "2_connected")
        return :random_2_connected
    end
    error("Invalid --infeas-base: $value. Expected globally_rigid, laman, or random_2_connected.")
end

function parse_contraction_vertices(value::AbstractString)
    normalized = strip(value)
    isempty(normalized) && return nothing, "auto"
    lowercase(normalized) == "auto" && return nothing, "auto"

    normalized = replace(normalized, "(" => "", ")" => "", " " => "")
    parts = split(normalized, ",")
    length(parts) == 2 || error("Invalid --contraction-vertices: $value. Expected auto or u,v.")
    u = parse_int(parts[1], "contraction vertex u")
    v = parse_int(parts[2], "contraction vertex v")
    return (u, v), "$u,$v"
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
    haskey(options, :type) || error("Missing required option --type")

    prefix, type_label = normalize_type(options[:type])
    target = abspath(get(options, :target, DEFAULT_TARGET))
    csv_path = haskey(options, :csv) ?
        abspath(options[:csv]) :
        joinpath(target, "instances.csv")
    contraction_vertices, contraction_vertices_label = parse_contraction_vertices(
        get(options, :contraction_vertices, "auto")
    )

    config = GeneratorConfig(
        prefix = prefix,
        type_label = type_label,
        count = parse_int(get(options, :count, "1"), "count"),
        target = target,
        csv_path = csv_path,
        instance_name_prefix = parse_instance_name_prefix(get(options, :instance_name_prefix, "")),
        n = parse_int(get(options, :n, "10"), "n"),
        R = parse_int(get(options, :R, "100"), "R"),
        num_anchors = parse_int(get(options, :num_anchors, "0"), "num_anchors"),
        alpha = parse_float(get(options, :alpha, "0.0"), "alpha"),
        seed_base = parse_int(get(options, :seed_base, "1"), "seed_base"),
        seed_step = parse_int(get(options, :seed_step, "1"), "seed_step"),
        edge_density = parse_float(get(options, :edge_density, "0.5"), "edge_density"),
        pH2 = parse_float(get(options, :pH2, "0.5"), "pH2"),
        max_coord_tries = parse_int(get(options, :max_coord_tries, "10000"), "max_coord_tries"),
        max_global_tries = parse_int(get(options, :max_global_tries, "10000"), "max_global_tries"),
        max_tries_H1 = parse_int(get(options, :max_tries_H1, "200"), "max_tries_H1"),
        max_tries_H2 = parse_int(get(options, :max_tries_H2, "300"), "max_tries_H2"),
        infeas_strategy = normalize_infeas_strategy(get(options, :infeas_strategy, "bounding_box")),
        infeas_base = normalize_infeas_base(get(options, :infeas_base, "globally_rigid")),
        box_margin = parse_int(get(options, :box_margin, "1"), "box_margin"),
        contraction_vertices = contraction_vertices,
        contraction_vertices_label = contraction_vertices_label,
    )

    config.count >= 1 || error("count must be >= 1")
    config.seed_step >= 1 || error("seed-step must be >= 1")
    config.box_margin >= 1 || error("box-margin must be >= 1")
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

function parse_instance_num(instance_name::AbstractString, instance_name_prefix::AbstractString, prefix::AbstractString)
    stem = "$(instance_name_prefix)$(prefix)_"
    startswith(instance_name, stem) || return nothing
    suffix = SubString(instance_name, ncodeunits(stem) + 1)
    return tryparse(Int, suffix)
end

function next_instance_number(
        csv_path::AbstractString,
        prefix::AbstractString,
        instance_name_prefix::AbstractString,
    )::Int
    !isfile(csv_path) && return 1
    filesize(csv_path) == 0 && return 1

    max_num = 0
    for (row_idx, row) in enumerate(CSV.File(csv_path))
        instance_name = row.instance_name === missing ? "" : string(row.instance_name)
        name_num = parse_instance_num(instance_name, instance_name_prefix, prefix)
        name_num === nothing && continue

        parsed_num = row.num === missing ? nothing : tryparse(Int, string(row.num))
        if parsed_num === nothing
            parsed_num = name_num
        end
        parsed_num === nothing && error(
            "Invalid num for prefix $(instance_name_prefix)$prefix on CSV row $(row_idx + 1) in $csv_path"
        )
        max_num = max(max_num, parsed_num)
    end

    return max_num + 1
end

function _same_bounding_box_family(row, config::GeneratorConfig)::Bool
    row.type === missing && return false
    row.n === missing && return false
    row.R === missing && return false
    row.num_anchors === missing && return false
    row.alpha === missing && return false
    row.infeas_strategy === missing && return false
    row.infeas_base === missing && return false
    row.box_margin === missing && return false
    common_matches = string(row.type) == config.type_label &&
        Int(row.n) == config.n &&
        Int(row.R) == config.R &&
        Int(row.num_anchors) == config.num_anchors &&
        Float64(row.alpha) == config.alpha &&
        string(row.infeas_strategy) == "bounding_box" &&
        string(row.infeas_base) == string(config.infeas_base) &&
        Int(row.box_margin) == config.box_margin
    common_matches || return false

    if config.infeas_base == :globally_rigid
        return row.max_tries_H2 !== missing && Int(row.max_tries_H2) == config.max_tries_H2
    elseif config.infeas_base == :laman
        return row.pH2 !== missing && Float64(row.pH2) == config.pH2 &&
            row.max_tries_H2 !== missing && Int(row.max_tries_H2) == config.max_tries_H2
    else
        return row.edge_density !== missing && Float64(row.edge_density) == config.edge_density &&
            row.max_coord_tries !== missing && Int(row.max_coord_tries) == config.max_coord_tries
    end
end

function previous_bounding_box_seed(config::GeneratorConfig)
    (!isfile(config.csv_path) || filesize(config.csv_path) == 0) && return nothing

    latest_seed = nothing
    for row in CSV.File(config.csv_path)
        instance_name = row.instance_name === missing ? "" : string(row.instance_name)
        parse_instance_num(instance_name, config.instance_name_prefix, config.prefix) === nothing &&
            continue
        _same_bounding_box_family(row, config) || continue
        row.seed === missing && error("Missing seed for existing bounding-box instance $instance_name")
        seed = Int(row.seed)
        latest_seed = latest_seed === nothing ? seed : max(latest_seed, seed)
    end
    return latest_seed
end

function base_kwargs(config::GeneratorConfig)
    if config.infeas_base == :globally_rigid
        return (
            R = config.R,
            max_global_tries = config.max_global_tries,
            max_tries_H2 = config.max_tries_H2,
        )
    elseif config.infeas_base == :laman
        return (
            R = config.R,
            pH2 = config.pH2,
            max_global_tries = config.max_global_tries,
            max_tries_H1 = config.max_tries_H1,
            max_tries_H2 = config.max_tries_H2,
        )
    elseif config.infeas_base == :random_2_connected
        return (
            R = config.R,
            edge_density = config.edge_density,
            max_coord_tries = config.max_coord_tries,
        )
    end
    error("Unsupported infeasible base: $(config.infeas_base)")
end

function generate_model(config::GeneratorConfig, seed::Int)
    if config.prefix == "con"
        return generate_2_connected_instance(
            config.n;
            R = config.R,
            edge_density = config.edge_density,
            seed = seed,
            max_coord_tries = config.max_coord_tries,
            num_anchors = config.num_anchors,
            alpha = config.alpha,
        )
    elseif config.prefix == "lam"
        return generate_laman_instance(
            config.n;
            R = config.R,
            pH2 = config.pH2,
            seed = seed,
            max_global_tries = config.max_global_tries,
            max_tries_H1 = config.max_tries_H1,
            max_tries_H2 = config.max_tries_H2,
            num_anchors = config.num_anchors,
            alpha = config.alpha,
        )
    elseif config.prefix == "gr"
        return generate_globally_rigid_instance(
            config.n;
            R = config.R,
            seed = seed,
            max_global_tries = config.max_global_tries,
            max_tries_H2 = config.max_tries_H2,
            num_anchors = config.num_anchors,
            alpha = config.alpha,
        )
    elseif config.prefix == "infeas"
        return generate_likely_infeasible_embedding_instance(
            config.n;
            strategy = config.infeas_strategy,
            base = config.infeas_base,
            seed = seed,
            num_anchors = config.num_anchors,
            alpha = config.alpha,
            box_margin = config.box_margin,
            contraction_vertices = config.contraction_vertices,
            base_kwargs(config)...,
        )
    end
    error("Unsupported type prefix: $(config.prefix)")
end

function generate_model_with_retry(config::GeneratorConfig, initial_seed::Int)
    candidate_seed = initial_seed
    for _ in 1:MAX_BOUNDING_BOX_CANDIDATE_TRIES
        try
            return generate_model(config, candidate_seed), candidate_seed
        catch err
            err isa BoundingBoxCandidateRejected || rethrow()
            println("rejected bounding-box seed $candidate_seed: $(sprint(showerror, err))")
            candidate_seed == typemax(Int) && error("cannot advance bounding-box candidate seed")
            candidate_seed += 1
        end
    end
    error(
        "failed to find a valid bounding-box candidate after " *
        "$MAX_BOUNDING_BOX_CANDIDATE_TRIES seeds starting at $initial_seed"
    )
end

function blank_csv_row()
    return Dict(column => "" for column in CSV_COLUMNS)
end

function set_common_row_fields!(
        row::Dict{String, String},
        config::GeneratorConfig,
        instance_name::String,
        num::Int,
        seed::Int,
    )
    row["instance_name"] = instance_name
    row["type"] = config.type_label
    row["num"] = string(num)
    row["created_at"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    row["n"] = string(config.n)
    row["R"] = string(config.R)
    row["seed"] = string(seed)
    row["num_anchors"] = string(config.num_anchors)
    row["alpha"] = string(config.alpha)
    return row
end

function set_type_specific_row_fields!(row::Dict{String, String}, config::GeneratorConfig)
    if config.prefix == "con"
        row["edge_density"] = string(config.edge_density)
        row["max_coord_tries"] = string(config.max_coord_tries)
    elseif config.prefix == "lam"
        row["pH2"] = string(config.pH2)
        row["max_tries_H2"] = string(config.max_tries_H2)
    elseif config.prefix == "gr"
        row["max_tries_H2"] = string(config.max_tries_H2)
    elseif config.prefix == "infeas"
        row["infeas_strategy"] = string(config.infeas_strategy)
        row["infeas_base"] = string(config.infeas_base)
        config.infeas_strategy == :bounding_box &&
            (row["box_margin"] = string(config.box_margin))
        if config.infeas_base == :globally_rigid
            row["max_tries_H2"] = string(config.max_tries_H2)
        elseif config.infeas_base == :laman
            row["pH2"] = string(config.pH2)
            row["max_tries_H2"] = string(config.max_tries_H2)
        elseif config.infeas_base == :random_2_connected
            row["edge_density"] = string(config.edge_density)
            row["max_coord_tries"] = string(config.max_coord_tries)
        end
    end
    return row
end

function ordered_csv_row(row::Dict{String, String})
    return (; (Symbol(column) => row[column] for column in CSV_COLUMNS)...)
end

function csv_row(config::GeneratorConfig, instance_name::String, num::Int, seed::Int)
    row = blank_csv_row()
    set_common_row_fields!(row, config, instance_name, num, seed)
    set_type_specific_row_fields!(row, config)
    return ordered_csv_row(row)
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

    start_num = next_instance_number(config.csv_path, config.prefix, config.instance_name_prefix)
    is_bounding_box = config.prefix == "infeas" && config.infeas_strategy == :bounding_box
    previous_seed = is_bounding_box ? previous_bounding_box_seed(config) : nothing
    next_bounding_box_seed = previous_seed === nothing ? nothing : previous_seed + config.seed_step
    for offset in 0:(config.count - 1)
        num = start_num + offset
        scheduled_seed = config.seed_base + (num - 1) * config.seed_step
        initial_seed = next_bounding_box_seed === nothing ?
            scheduled_seed : max(scheduled_seed, next_bounding_box_seed)
        instance_name = "$(config.instance_name_prefix)$(config.prefix)_$num"
        file_name = "$instance_name.lp"
        file_path = joinpath(config.target, file_name)
        isfile(file_path) && error("Refusing to overwrite existing instance file: $file_path")

        generated, seed = is_bounding_box ?
            generate_model_with_retry(config, initial_seed) :
            (generate_model(config, scheduled_seed), scheduled_seed)
        model, _, _ = generated
        is_bounding_box && (next_bounding_box_seed = seed + config.seed_step)
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
