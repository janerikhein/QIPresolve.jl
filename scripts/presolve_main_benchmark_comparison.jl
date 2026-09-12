#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
end

module MainBenchmarkPresolveComparison

using CSV
using Statistics: mean
import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC
import SCIP

include("residue_constraint_experiment_common.jl")
const Common = ResidueConstraintExperimentCommon
const ROOT = normpath(joinpath(@__DIR__, ".."))
const METRICS = (:domain_reduction_iso, :bound_tightening_iso,
    :domain_reduction_comb, :bound_tightening_comb, :domain_reduction_scip,
    :bound_tightening_scip, :domain_reduction_full, :bound_tightening_full)

Base.@kwdef struct Config
    input_dir::String = joinpath(ROOT, "instances", "main_benchmark")
    output_dir::String = joinpath(ROOT, "results", "main_benchmark_presolve_comparison")
    parity_strategy::Symbol = QIP.PresolveConfig.DEFAULT_PRESOLVE_PARITY_STRATEGY
    residue_strategy::Symbol = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY
    residue_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD
    treewidth_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD
    scip_config::Union{Nothing, String} = nothing
    limit::Int = typemax(Int)
end

function usage()
    return """
    Usage: julia --project=. scripts/presolve_main_benchmark_comparison.jl [options]
      --input-dir DIR             Default: instances/main_benchmark
      --output-dir DIR            Default: results/main_benchmark_presolve_comparison
      --parity-strategy NAME      full, mod2_basic, mod4_basic
      --residue-strategy NAME     divisor_free, small_primes, powers_of_two
      --residue-threshold N       Default: $(QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD)
      --treewidth-threshold N     Default: $(QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD)
      --scip-config FILE          Optional .set file for both SCIP runs
      --limit N                  Process the first N instances (smoke runs)
      -h, --help                 Show this help

    SCIP uses presolving only, with default settings unless overridden.
    Removed original constraints contribute 1.0 to bound tightening.
    Outputs: per_instance.csv, aggregated.csv, diagnostics.log, run_config.txt
    """
end

function parse_args(args::Vector{String})
    options = Dict{Symbol, Any}()
    i = 1
    while i <= length(args)
        arg = args[i]
        arg in ("-h", "--help") && (println(usage()); return nothing)
        startswith(arg, "--") || error("Unexpected argument: $arg")
        parts = split(arg[3:end], '='; limit = 2)
        if length(parts) == 1
            i += 1
            i <= length(args) || error("Missing value for $arg")
            push!(parts, args[i])
        end
        key = Symbol(replace(parts[1], '-' => '_'))
        key in fieldnames(Config) || error("Unknown option: $(parts[1])")
        value = String(parts[2])
        options[key] = if key in (:limit, :residue_threshold, :treewidth_threshold)
            parse(Int, value)
        elseif key in (:parity_strategy, :residue_strategy)
            Symbol(replace(value, '-' => '_'))
        else
            abspath(value)
        end
        i += 1
    end
    config = Config(; options...)
    config.limit > 0 || error("--limit must be positive")
    config.residue_threshold >= 0 || error("--residue-threshold must be nonnegative")
    config.treewidth_threshold >= 0 || error("--treewidth-threshold must be nonnegative")
    PC._normalize_parity_strategy(config.parity_strategy)
    PC._generate_residue_moduli(config.residue_strategy, config.residue_threshold)
    return config
end

"Reassemble complementary LP rows without performing variable or GCD reductions."
function rejoin_ranges!(model::PC.QPModel)
    lower_rows = Dict{Any, Vector{PC.Constraint}}()
    upper_rows = Dict{Any, Vector{PC.Constraint}}()
    removed = Base.IdSet{PC.Constraint}()
    for con in model.cons
        key = PC._constraint_coefficient_key(con)
        # Canonicalize the sign so f <= b and -f <= -a also pair.
        coeffs = [last(t) for t in key.lin_terms]
        append!(coeffs, [last(t) for t in key.quad_terms])
        isempty(coeffs) && continue
        negative = first(coeffs) < 0
        canonical = negative ? PC._negated_constraint_key(key) : key
        lhs, rhs = negative ? (-con.rhs, -con.lhs) : (con.lhs, con.rhs)
        if isfinite(lhs) && rhs == Inf
            push!(get!(lower_rows, canonical, PC.Constraint[]), con)
        elseif lhs == -Inf && isfinite(rhs)
            push!(get!(upper_rows, canonical, PC.Constraint[]), con)
        end
    end
    for (key, lowers) in lower_rows
        uppers = get(upper_rows, key, PC.Constraint[])
        for (lower, upper) in zip(lowers, uppers)
            same_sign = PC._constraint_coefficient_key(lower) == PC._constraint_coefficient_key(upper)
            lhs, rhs = same_sign ? (upper.lhs, upper.rhs) : (-upper.rhs, -upper.lhs)
            lower.lhs = max(lower.lhs, lhs)
            lower.rhs = min(lower.rhs, rhs)
            push!(removed, upper)
        end
    end
    filter!(c -> !(c in removed), model.cons)
    return model
end

function load_instance(path::AbstractString)
    lp = FF.Model(format = FF.FORMAT_LP)
    MOI.read_from_file(lp, path)
    model = QIP.build_model(QIP.from_moi(lp))
    rejoin_ranges!(model)
    model.infeasible |= any(c -> c.lhs > c.rhs, model.cons)
    for var in values(model.vars)
        isfinite(var.lb) && isfinite(var.ub) || error("Domain metric requires finite variable bounds: $path")
    end
    return model
end

log_domain_sum(model::PC.QPModel) = sum(
    log(v.ub - v.lb + 1.0) for v in values(model.vars); init = 0.0)
domain_reduction(before::Real, after::Real) = before == 0.0 ? 0.0 : (before - after) / before

struct BoundBaseline
    con::PC.Constraint
    lhs::Float64
    rhs::Float64
    scale::Float64
end
bound_baselines(model) = [BoundBaseline(c, c.lhs, c.rhs, c._bound_scale) for c in model.cons]

function bound_contribution(baseline::BoundBaseline, after_width::Real; removed::Bool = false)
    removed && return 1.0
    # Invoke the reference helper with translated bounds: only width matters.
    before = (lhs = baseline.lhs, rhs = baseline.rhs)
    con = PC.Constraint(0, PC.QuadExpr(Tuple{Float64,Int,Int}[], Tuple{Float64,Int}[]), 0.0, 0.0)
    con.lhs = 0.0
    con.rhs = Float64(after_width)
    return Common.relative_bound_range_reduction(before, con)
end

function core_contributions(model, baselines)
    alive = union!(Base.IdSet{PC.Constraint}(), model.cons)
    return [bound_contribution(b, (b.con.rhs - b.con.lhs) * b.scale / b.con._bound_scale;
        removed = !(b.con in alive)) for b in baselines]
end

average_contributions(values) = isempty(values) ? 0.0 : mean(values)
core_status(model) = model.infeasible ? "infeasible" : isempty(model.cons) ? "feasible" : "reduced"

function run_core(baseline, config; enable_parity, enable_residue)
    model = deepcopy(baseline)
    bounds = bound_baselines(model)
    result = QIP.presolve!(model; enable_parity, enable_residue,
        parity_strategy = config.parity_strategy,
        residue_strategy = config.residue_strategy,
        residue_threshold = config.residue_threshold,
        treewidth_threshold = config.treewidth_threshold)
    return (result = result, bounds = bounds, contributions = core_contributions(model, bounds),
        status = core_status(model), log_domain = model.infeasible ? 0.0 : log_domain_sum(model))
end

include("presolve_comparison_scip.jl")

function with_strategy(f::Function, instance, strategy)
    try
        return f()
    catch err
        error("Instance $instance, strategy $strategy: $(sprint(showerror, err))")
    end
end

function evaluate_instance(path, type, config; diagnostics::IO = stderr)
    instance = splitext(basename(path))[1]
    original = with_strategy(instance, "load") do
        load_instance(path)
    end
    original_log = log_domain_sum(original)
    parity = with_strategy(instance, "parity") do
        run_core(original, config; enable_parity = true, enable_residue = false)
    end
    residue = with_strategy(instance, "residue") do
        run_core(original, config; enable_parity = false, enable_residue = true)
    end
    combined = with_strategy(instance, "combined") do
        run_core(original, config; enable_parity = true, enable_residue = true)
    end
    scip = with_strategy(instance, "SCIP") do
        run_scip(original, bound_baselines(original), config; diagnostics, label = "$instance/scip")
    end
    full = if combined.status != "reduced"
        combined
    else
        with_strategy(instance, "full") do
            run_scip(combined.result.model, combined.bounds, config;
                diagnostics, label = "$instance/full")
        end
    end
    return (instance_name = instance, type = type,
        domain_reduction_iso = domain_reduction(original_log, parity.log_domain),
        bound_tightening_iso = average_contributions(residue.contributions),
        domain_reduction_comb = domain_reduction(original_log, combined.log_domain),
        bound_tightening_comb = average_contributions(combined.contributions),
        domain_reduction_scip = domain_reduction(original_log, scip.log_domain),
        bound_tightening_scip = average_contributions(scip.contributions),
        domain_reduction_full = domain_reduction(original_log, full.log_domain),
        bound_tightening_full = average_contributions(full.contributions),
        status_parity = parity.status, status_residue = residue.status,
        status_comb = combined.status, status_scip = scip.status, status_full = full.status)
end

function instance_types(directory)
    types = Dict{String, String}()
    for filename in ("random_instances.csv", "embedding_instances.csv")
        path = joinpath(directory, filename)
        isfile(path) || continue
        for row in CSV.File(path)
            file = String(row.file_name)
            haskey(types, file) && error("Duplicate metadata for $file")
            types[file] = if filename == "random_instances.csv"
                "random_$(row.subtype)"
            elseif hasproperty(row, :infeas_strategy) && !ismissing(row.infeas_strategy)
                "embedding_$(row.exactness)_infeasible_$(row.infeas_strategy)"
            else
                "embedding_$(row.exactness)_$(row.anchoring)_$(row.graph_type)"
            end
        end
    end
    return types
end

function aggregate_rows(rows)
    return [merge((type = type, instance_count = count(r -> r.type == type, rows)),
        NamedTuple{METRICS}(Tuple(mean(getproperty(r, key) for r in rows if r.type == type)
            for key in METRICS))) for type in sort!(unique([r.type for r in rows]))]
end

function run_experiment(config::Config)
    isdir(config.input_dir) || error("Input directory not found: $(config.input_dir)")
    config.scip_config === nothing || isfile(config.scip_config) || error("SCIP settings file not found")
    files = sort!(filter(f -> endswith(lowercase(f), ".lp"), readdir(config.input_dir)))
    isempty(files) && error("No LP instances found in $(config.input_dir)")
    resize!(files, min(length(files), config.limit))
    types = instance_types(config.input_dir)
    all(f -> haskey(types, f), files) || error("Missing benchmark metadata for one or more LP files")
    mkpath(config.output_dir)
    open(joinpath(config.output_dir, "run_config.txt"), "w") do io
        println(io, "Julia: ", VERSION)
        println(io, "SCIP: ", SCIP.SCIPmajorVersion(), '.', SCIP.SCIPminorVersion(), '.', SCIP.SCIPtechVersion())
        println(io, "Instances: ", length(files))
        for key in fieldnames(Config)
            println(io, key, " = ", getfield(config, key))
        end
        config.scip_config === nothing || print(io, "\nSCIP settings contents:\n", read(config.scip_config, String))
    end
    rows = NamedTuple[]
    open(joinpath(config.output_dir, "diagnostics.log"), "w") do diagnostics
        for (index, file) in enumerate(files)
            println("[$index/$(length(files))] $file")
            flush(stdout)
            row = evaluate_instance(joinpath(config.input_dir, file), types[file], config; diagnostics)
            push!(rows, row)
            CSV.write(joinpath(config.output_dir, "per_instance.csv"), [row]; append = index > 1)
            CSV.write(joinpath(config.output_dir, "aggregated.csv"), aggregate_rows(rows))
            flush(diagnostics)
        end
    end
    println("Wrote $(length(rows)) instances and $(length(aggregate_rows(rows))) types to $(config.output_dir)")
    return rows
end

function main(args = copy(ARGS))
    config = parse_args(args)
    config === nothing || run_experiment(config)
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    MainBenchmarkPresolveComparison.main()
end
