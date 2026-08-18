#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module PresolveLpScipStatsScript

using Dates
import JSON
import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC
import SCIP

const DEFAULT_VALIDATION_TOLERANCE = 1.0e-6

const SCIP_METRIC_KEYS = [
    "label",
    "skipped",
    "error",
    "solver_name",
    "solver_version",
    "termination_status",
    "primal_status",
    "dual_status",
    "raw_status",
    "result_count",
    "objective_value",
    "objective_bound",
    "relative_gap",
    "solve_time_sec",
    "wall_time_sec",
    "node_count",
    "simplex_iterations",
    "n_solutions",
    "n_solutions_found",
    "n_limit_solutions_found",
    "n_best_solutions_found",
    "n_original_variables",
    "n_original_binary_variables",
    "n_original_integer_variables",
    "n_original_implicit_integer_variables",
    "n_original_continuous_variables",
    "n_original_constraints",
    "n_transformed_variables",
    "n_transformed_binary_variables",
    "n_transformed_integer_variables",
    "n_transformed_implicit_integer_variables",
    "n_transformed_continuous_variables",
    "n_transformed_constraints",
    "n_active_constraints",
    "n_enabled_constraints",
    "n_lps",
    "n_lp_iterations",
    "n_root_lp_iterations",
    "n_node_lp_iterations",
    "n_cuts",
    "n_cuts_found",
    "n_cuts_applied",
    "n_conflict_constraints_found",
    "n_conflict_constraints_applied",
    "metric_errors",
]

Base.@kwdef struct CliConfig
    input_path::String
    scip_config::Union{Nothing, String} = nothing
    output_path::Union{Nothing, String} = nothing
    out_dir::Union{Nothing, String} = nothing
    time_limit_sec::Union{Nothing, Float64} = nothing
    residue_strategy::Symbol = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY
    residue_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD
    treewidth_threshold::Int = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD
    silent::Bool = false
    validation_tolerance::Float64 = DEFAULT_VALIDATION_TOLERANCE
end

function usage()
    return """
    Usage:
      julia --project=. scripts/presolve_lp_scip_stats.jl instance.lp [options]

    Options:
      --scip-config file.set       SCIP parameter file applied to both SCIP runs
      --output path                Explicit JSON output path
      --out-dir dir                Directory for the default <instance>_stats.json
      --time-limit sec             SCIP time limit in seconds
      --residue-strategy name      small_primes, powers_of_two, or divisor_free
      --residue-threshold n        Maximum residue modulus threshold
      --treewidth-threshold n      Residue DP treewidth threshold
      --silent                     Set SCIP display/verblevel to 0 after config
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

function _parse_float(value::AbstractString, name::AbstractString)
    try
        parsed = parse(Float64, value)
        parsed >= 0.0 || error()
        return parsed
    catch
        error("Invalid nonnegative Float64 for $name: $value")
    end
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
    scip_config = nothing
    output_path = nothing
    out_dir = nothing
    time_limit_sec = nothing
    residue_strategy = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY
    residue_threshold = QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD
    treewidth_threshold = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD
    silent = false

    index = 1
    while index <= length(args)
        arg = args[index]
        if arg == "--silent"
            silent = true
            index += 1
        elseif startswith(arg, "--")
            option, value, consumed = _option_value(args, index, arg)
            if option == "--scip-config"
                scip_config = value
            elseif option == "--output"
                output_path = value
            elseif option == "--out-dir"
                out_dir = value
            elseif option == "--time-limit"
                time_limit_sec = _parse_float(value, option)
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
        scip_config = scip_config === nothing ? nothing : abspath(scip_config),
        output_path = output_path === nothing ? nothing : abspath(output_path),
        out_dir = out_dir === nothing ? nothing : abspath(out_dir),
        time_limit_sec = time_limit_sec,
        residue_strategy = residue_strategy,
        residue_threshold = residue_threshold,
        treewidth_threshold = treewidth_threshold,
        silent = silent,
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

function _has_quadratic_terms(expr::PC.QuadExpr)
    var_ids = collect(PC.vars(expr))
    for (index, var_id_1) in enumerate(var_ids)
        for var_id_2 in @view var_ids[index:end]
            isapprox(PC.get_quad_coeff(expr, var_id_1, var_id_2), 0.0) || return true
        end
    end
    return false
end

function ensure_scip_supported!(model::PC.QPModel)
    if _has_quadratic_terms(model.obj_expr)
        throw(
            ArgumentError(
                "Unsupported LP: the SCIP-backed MOI path currently supports affine objectives only; " *
                "quadratic objectives cannot be solved by this script.",
            ),
        )
    end
    return model
end

function scip_optimizer(moi_model)
    moi_model isa SCIP.Optimizer && return moi_model
    if hasfield(typeof(moi_model), :model)
        optimizer = getfield(moi_model, :model)
        optimizer isa SCIP.Optimizer && return optimizer
    end
    error("Unable to locate SCIP.Optimizer backend in $(typeof(moi_model))")
end

function apply_scip_settings!(optimizer::SCIP.Optimizer, config::CliConfig)
    if config.scip_config !== nothing
        isfile(config.scip_config) || error("SCIP config file not found: $(config.scip_config)")
        SCIP.@SCIP_CALL SCIP.SCIPreadParams(optimizer.inner, config.scip_config)
    end

    if config.time_limit_sec !== nothing
        MOI.set(optimizer, MOI.TimeLimitSec(), config.time_limit_sec)
    end

    config.silent && MOI.set(optimizer, MOI.Silent(), true)
    return optimizer
end

function build_configured_scip_model(model::PC.QPModel, config::CliConfig)
    moi_model = QIP.build_moi_model(model, :scip)
    apply_scip_settings!(scip_optimizer(moi_model), config)
    return moi_model
end

function empty_scip_metrics(label::AbstractString; skipped::Bool = false, error = nothing)
    metrics = Dict{String, Any}(key => nothing for key in SCIP_METRIC_KEYS)
    metrics["label"] = label
    metrics["skipped"] = skipped
    metrics["error"] = error
    metrics["metric_errors"] = Dict{String, Any}()
    return metrics
end

function skipped_scip_metrics(label::AbstractString, reason::AbstractString)
    return empty_scip_metrics(label; skipped = true, error = reason)
end

function _record_metric!(metrics::Dict{String, Any}, key::AbstractString, f::Function)
    errors = metrics["metric_errors"]::Dict{String, Any}
    try
        metrics[key] = json_safe(f())
    catch err
        metrics[key] = nothing
        errors[key] = sprint(showerror, err)
    end
    return metrics
end

record_metric!(metrics::Dict{String, Any}, key::AbstractString, f::Function) =
    _record_metric!(metrics, key, f)

record_metric!(f::Function, metrics::Dict{String, Any}, key::AbstractString) =
    _record_metric!(metrics, key, f)

function record_string_metric!(metrics::Dict{String, Any}, key::AbstractString, f::Function)
    return _record_metric!(metrics, key, () -> string(f()))
end

record_string_metric!(f::Function, metrics::Dict{String, Any}, key::AbstractString) =
    record_string_metric!(metrics, key, f)

function _safe_result_count(moi_model)
    try
        return MOI.get(moi_model, MOI.ResultCount())
    catch
        return 0
    end
end

function collect_scip_metrics(moi_model, label::AbstractString, wall_time_sec::Float64)
    metrics = empty_scip_metrics(label)
    optimizer = scip_optimizer(moi_model)
    scip = optimizer.inner

    metrics["wall_time_sec"] = wall_time_sec

    record_string_metric!(metrics, "solver_name") do
        MOI.get(moi_model, MOI.SolverName())
    end
    record_string_metric!(metrics, "solver_version") do
        MOI.get(moi_model, MOI.SolverVersion())
    end
    record_string_metric!(metrics, "termination_status") do
        MOI.get(moi_model, MOI.TerminationStatus())
    end
    record_string_metric!(metrics, "raw_status") do
        MOI.get(moi_model, MOI.RawStatusString())
    end
    record_metric!(metrics, "result_count") do
        MOI.get(moi_model, MOI.ResultCount())
    end

    result_count = metrics["result_count"] isa Integer ? metrics["result_count"] : _safe_result_count(moi_model)
    if result_count > 0
        record_string_metric!(metrics, "primal_status") do
            MOI.get(moi_model, MOI.PrimalStatus())
        end
        record_metric!(metrics, "objective_value") do
            MOI.get(moi_model, MOI.ObjectiveValue())
        end
    else
        metrics["primal_status"] = string(MOI.NO_SOLUTION)
        metrics["objective_value"] = nothing
    end

    record_string_metric!(metrics, "dual_status") do
        MOI.get(moi_model, MOI.DualStatus())
    end
    record_metric!(metrics, "objective_bound") do
        MOI.get(moi_model, MOI.ObjectiveBound())
    end
    record_metric!(metrics, "relative_gap") do
        MOI.get(moi_model, MOI.RelativeGap())
    end
    record_metric!(metrics, "solve_time_sec") do
        MOI.get(moi_model, MOI.SolveTimeSec())
    end
    record_metric!(metrics, "node_count") do
        MOI.get(moi_model, MOI.NodeCount())
    end
    record_metric!(metrics, "simplex_iterations") do
        MOI.get(moi_model, MOI.SimplexIterations())
    end

    record_metric!(metrics, "n_solutions") do
        SCIP.SCIPgetNSols(scip)
    end
    record_metric!(metrics, "n_solutions_found") do
        SCIP.SCIPgetNSolsFound(scip)
    end
    record_metric!(metrics, "n_limit_solutions_found") do
        SCIP.SCIPgetNLimSolsFound(scip)
    end
    record_metric!(metrics, "n_best_solutions_found") do
        SCIP.SCIPgetNBestSolsFound(scip)
    end

    record_metric!(metrics, "n_original_variables") do
        SCIP.SCIPgetNOrigVars(scip)
    end
    record_metric!(metrics, "n_original_binary_variables") do
        SCIP.SCIPgetNOrigBinVars(scip)
    end
    record_metric!(metrics, "n_original_integer_variables") do
        SCIP.SCIPgetNOrigIntVars(scip)
    end
    record_metric!(metrics, "n_original_implicit_integer_variables") do
        SCIP.SCIPgetNOrigImplVars(scip)
    end
    record_metric!(metrics, "n_original_continuous_variables") do
        SCIP.SCIPgetNOrigContVars(scip)
    end
    record_metric!(metrics, "n_original_constraints") do
        SCIP.SCIPgetNOrigConss(scip)
    end

    record_metric!(metrics, "n_transformed_variables") do
        SCIP.SCIPgetNVars(scip)
    end
    record_metric!(metrics, "n_transformed_binary_variables") do
        SCIP.SCIPgetNBinVars(scip)
    end
    record_metric!(metrics, "n_transformed_integer_variables") do
        SCIP.SCIPgetNIntVars(scip)
    end
    record_metric!(metrics, "n_transformed_implicit_integer_variables") do
        SCIP.SCIPgetNImplVars(scip)
    end
    record_metric!(metrics, "n_transformed_continuous_variables") do
        SCIP.SCIPgetNContVars(scip)
    end
    record_metric!(metrics, "n_transformed_constraints") do
        SCIP.SCIPgetNConss(scip)
    end
    record_metric!(metrics, "n_active_constraints") do
        SCIP.SCIPgetNActiveConss(scip)
    end
    record_metric!(metrics, "n_enabled_constraints") do
        SCIP.SCIPgetNEnabledConss(scip)
    end

    record_metric!(metrics, "n_lps") do
        SCIP.SCIPgetNLPs(scip)
    end
    record_metric!(metrics, "n_lp_iterations") do
        SCIP.SCIPgetNLPIterations(scip)
    end
    record_metric!(metrics, "n_root_lp_iterations") do
        SCIP.SCIPgetNRootLPIterations(scip)
    end
    record_metric!(metrics, "n_node_lp_iterations") do
        SCIP.SCIPgetNNodeLPIterations(scip)
    end

    record_metric!(metrics, "n_cuts") do
        SCIP.SCIPgetNCuts(scip)
    end
    record_metric!(metrics, "n_cuts_found") do
        SCIP.SCIPgetNCutsFound(scip)
    end
    record_metric!(metrics, "n_cuts_applied") do
        SCIP.SCIPgetNCutsApplied(scip)
    end
    record_metric!(metrics, "n_conflict_constraints_found") do
        SCIP.SCIPgetNConflictConssFound(scip)
    end
    record_metric!(metrics, "n_conflict_constraints_applied") do
        SCIP.SCIPgetNConflictConssApplied(scip)
    end

    return metrics
end

function variable_id_from_name(name::AbstractString)
    matched = match(r"^x([0-9]+)$", name)
    matched === nothing && return nothing
    return parse(Int, matched.captures[1])
end

function reduced_solution(moi_model)
    solution = Dict{PC.VarId, Float64}()
    variables = MOI.get(moi_model, MOI.ListOfVariableIndices())

    for variable in variables
        name = MOI.get(moi_model, MOI.VariableName(), variable)
        var_id = variable_id_from_name(name)
        var_id === nothing && continue
        solution[var_id] = Float64(MOI.get(moi_model, MOI.VariablePrimal(), variable))
    end

    return solution
end

function _feasible_primal_status(metrics::Dict{String, Any})
    status = get(metrics, "primal_status", nothing)
    return status == string(MOI.FEASIBLE_POINT) || status == string(MOI.NEARLY_FEASIBLE_POINT)
end

function _metric_int(metrics::Dict{String, Any}, key::AbstractString, default::Int = 0)
    value = get(metrics, key, default)
    return value isa Integer ? Int(value) : default
end

function solve_with_scip(model::PC.QPModel, label::AbstractString, config::CliConfig)
    moi_model = nothing
    try
        moi_model = build_configured_scip_model(model, config)
    catch err
        return (
            metrics = empty_scip_metrics(label; error = sprint(showerror, err)),
            reduced_solution = Dict{PC.VarId, Float64}(),
        )
    end

    start_time = time_ns()
    optimize_error = nothing
    try
        MOI.optimize!(moi_model)
    catch err
        optimize_error = sprint(showerror, err)
    end
    wall_time_sec = (time_ns() - start_time) / 1.0e9

    metrics = collect_scip_metrics(moi_model, label, wall_time_sec)
    metrics["error"] = optimize_error

    solution = Dict{PC.VarId, Float64}()
    if optimize_error === nothing && _metric_int(metrics, "result_count") > 0 && _feasible_primal_status(metrics)
        try
            solution = reduced_solution(moi_model)
        catch err
            metrics["metric_errors"]["reduced_solution"] = sprint(showerror, err)
        end
    end

    return (metrics = metrics, reduced_solution = solution)
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

function evaluate_expr(expr::PC.QuadExpr, solution::Dict{PC.VarId, Float64})
    total = expr.constant
    var_ids = collect(PC.vars(expr))

    for var_id in var_ids
        haskey(solution, var_id) || error("missing solution value for variable $var_id")
        total += PC.get_lin_coeff(expr, var_id) * solution[var_id]
    end

    for (index, var_id_1) in enumerate(var_ids)
        x1 = solution[var_id_1]
        for var_id_2 in @view var_ids[index:end]
            coeff = PC.get_quad_coeff(expr, var_id_1, var_id_2)
            coeff == 0.0 && continue
            total += coeff * x1 * solution[var_id_2]
        end
    end

    return Float64(total)
end

function _bound_violation(value::Float64, lb::Float64, ub::Float64)
    lower = isfinite(lb) ? max(0.0, lb - value) : 0.0
    upper = isfinite(ub) ? max(0.0, value - ub) : 0.0
    return max(lower, upper)
end

function _constraint_violation(activity::Float64, lhs::Float64, rhs::Float64)
    lower = isfinite(lhs) ? max(0.0, lhs - activity) : 0.0
    upper = isfinite(rhs) ? max(0.0, activity - rhs) : 0.0
    return max(lower, upper)
end

function solution_dict(solution::Dict{PC.VarId, Float64})
    result = Dict{String, Any}()
    for var_id in sort!(collect(keys(solution)))
        result["x$(var_id)"] = json_safe(solution[var_id])
    end
    return result
end

function skipped_validation(reason::AbstractString; tolerance::Float64 = DEFAULT_VALIDATION_TOLERANCE)
    return Dict{String, Any}(
        "checked" => false,
        "feasible" => false,
        "status" => "skipped",
        "reason" => reason,
        "tolerance" => tolerance,
        "reduced_solution_size" => 0,
        "postsolved_solution_size" => 0,
        "missing_original_values" => 0,
        "variable_bound_violation_count" => 0,
        "variable_integrality_violation_count" => 0,
        "constraint_violation_count" => 0,
        "max_variable_bound_violation" => 0.0,
        "max_variable_integrality_violation" => 0.0,
        "max_constraint_violation" => 0.0,
        "max_violation" => 0.0,
        "original_objective_sense" => nothing,
        "original_objective_value" => nothing,
        "reduced_solution" => Dict{String, Any}(),
        "postsolved_solution" => Dict{String, Any}(),
        "errors" => String[],
    )
end

function validate_postsolved_solution(
        original_model::PC.QPModel,
        postsolver::PC.ParityPostsolver,
        reduced::Dict{PC.VarId, Float64},
        solve_metrics::Dict{String, Any};
        tolerance::Float64,
    )
    if _metric_int(solve_metrics, "result_count") <= 0 || !_feasible_primal_status(solve_metrics)
        return skipped_validation(
            "Presolved SCIP run did not return a feasible primal solution";
            tolerance = tolerance,
        )
    end

    errors = String[]
    original = Dict{PC.VarId, Float64}()
    try
        original = PC.postsolve(postsolver, reduced)
    catch err
        push!(errors, sprint(showerror, err))
        validation = skipped_validation("Postsolve failed"; tolerance = tolerance)
        validation["checked"] = true
        validation["status"] = "postsolve_error"
        validation["reduced_solution_size"] = length(reduced)
        validation["reduced_solution"] = solution_dict(reduced)
        validation["errors"] = errors
        return validation
    end

    missing_original_values = 0
    variable_bound_violation_count = 0
    variable_integrality_violation_count = 0
    constraint_violation_count = 0
    max_variable_bound_violation = 0.0
    max_variable_integrality_violation = 0.0
    max_constraint_violation = 0.0

    for (var_id, var) in original_model.vars
        if !haskey(original, var_id)
            missing_original_values += 1
            continue
        end

        value = original[var_id]
        bound_violation = _bound_violation(value, var.lb, var.ub)
        max_variable_bound_violation = max(max_variable_bound_violation, bound_violation)
        bound_violation > tolerance && (variable_bound_violation_count += 1)

        integrality_violation = abs(value - round(value))
        max_variable_integrality_violation = max(max_variable_integrality_violation, integrality_violation)
        integrality_violation > tolerance && (variable_integrality_violation_count += 1)
    end

    for con in original_model.cons
        try
            activity = evaluate_expr(con.qe, original)
            violation = _constraint_violation(activity, con.lhs, con.rhs)
            max_constraint_violation = max(max_constraint_violation, violation)
            violation > tolerance && (constraint_violation_count += 1)
        catch err
            push!(errors, sprint(showerror, err))
            constraint_violation_count += 1
            max_constraint_violation = Inf
        end
    end

    objective_value = nothing
    try
        objective_value = evaluate_expr(original_model.obj_expr, original)
    catch err
        push!(errors, sprint(showerror, err))
    end

    max_violation = max(
        max_variable_bound_violation,
        max_variable_integrality_violation,
        max_constraint_violation,
    )
    feasible =
        missing_original_values == 0 &&
        variable_bound_violation_count == 0 &&
        variable_integrality_violation_count == 0 &&
        constraint_violation_count == 0 &&
        isempty(errors) &&
        max_violation <= tolerance

    return Dict{String, Any}(
        "checked" => true,
        "feasible" => feasible,
        "status" => feasible ? "feasible" : "violated",
        "reason" => nothing,
        "tolerance" => tolerance,
        "reduced_solution_size" => length(reduced),
        "postsolved_solution_size" => length(original),
        "missing_original_values" => missing_original_values,
        "variable_bound_violation_count" => variable_bound_violation_count,
        "variable_integrality_violation_count" => variable_integrality_violation_count,
        "constraint_violation_count" => constraint_violation_count,
        "max_variable_bound_violation" => json_safe(max_variable_bound_violation),
        "max_variable_integrality_violation" => json_safe(max_variable_integrality_violation),
        "max_constraint_violation" => json_safe(max_constraint_violation),
        "max_violation" => json_safe(max_violation),
        "original_objective_sense" => string(original_model.obj_sense),
        "original_objective_value" => json_safe(objective_value),
        "reduced_solution" => solution_dict(reduced),
        "postsolved_solution" => solution_dict(original),
        "errors" => errors,
    )
end

function scip_version_string()
    try
        return string(
            SCIP.SCIPmajorVersion(),
            ".",
            SCIP.SCIPminorVersion(),
            ".",
            SCIP.SCIPtechVersion(),
        )
    catch err
        return sprint(showerror, err)
    end
end

function parameters_dict(config::CliConfig)
    return Dict{String, Any}(
        "scip_config" => config.scip_config,
        "time_limit_sec" => config.time_limit_sec,
        "residue_strategy" => string(config.residue_strategy),
        "residue_threshold" => config.residue_threshold,
        "treewidth_threshold" => config.treewidth_threshold,
        "silent" => config.silent,
        "validation_tolerance" => config.validation_tolerance,
    )
end

function metadata_dict(config::CliConfig, output_path::AbstractString, load_info)
    return Dict{String, Any}(
        "generated_at_utc" => _utc_timestamp(),
        "input_path" => config.input_path,
        "input_size_bytes" => filesize(config.input_path),
        "output_path" => output_path,
        "julia_version" => string(VERSION),
        "scip_version" => scip_version_string(),
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
    config.scip_config === nothing || isfile(config.scip_config) ||
        error("SCIP config file not found: $(config.scip_config)")

    output_path = default_output_path(config)
    instance_name = splitext(basename(config.input_path))[1]

    load_info = load_lp_core_model(config.input_path)
    loaded_model = load_info.model
    ensure_scip_supported!(loaded_model)

    model_loaded = PC.ModelStats(loaded_model)
    original_solve = solve_with_scip(deepcopy(loaded_model), "original", config)

    presolve_model = deepcopy(loaded_model)
    presolve_result = run_tracked_presolve!(presolve_model, config)

    if presolve_result.model.infeasible
        presolved_metrics = skipped_scip_metrics("presolved", "QIPresolve presolve proved the model infeasible")
        presolved_solution = Dict{PC.VarId, Float64}()
        validation = skipped_validation(
            "QIPresolve presolve proved the model infeasible";
            tolerance = config.validation_tolerance,
        )
    else
        presolved_solve = solve_with_scip(deepcopy(presolve_result.model), "presolved", config)
        presolved_metrics = presolved_solve.metrics
        presolved_solution = presolved_solve.reduced_solution
        validation = validate_postsolved_solution(
            loaded_model,
            presolve_result.postsolver,
            presolved_solution,
            presolved_metrics;
            tolerance = config.validation_tolerance,
        )
    end

    data = Dict{String, Any}(
        "instance_name" => instance_name,
        "metadata" => metadata_dict(config, output_path, load_info),
        "parameters" => parameters_dict(config),
        "scip_original" => original_solve.metrics,
        "scip_presolved" => presolved_metrics,
        "model_loaded" => struct_stats_dict(model_loaded),
        "model_before" => struct_stats_dict(presolve_result.model_before),
        "model_after" => struct_stats_dict(presolve_result.model_after),
        "parity_stats" => struct_stats_dict(presolve_result.parity_stats),
        "residue_stats" => struct_stats_dict(presolve_result.residue_stats),
        "validation" => validation,
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
    PresolveLpScipStatsScript.main(copy(ARGS))
end
