#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Graphs
using JuMP: backend
import MathOptInterface as MOI
import QIPresolve as QIP

PC = QIP.PresolvingCore
GE = QIP.InstanceGeneration

const N = 10
const R = 50
const PH2 = 0.5
const SEED = 17
const USE_PARITY_PRESOLVE = false
const TIME_LIMIT_SEC = nothing
const SCIP_VERBLEVEL = 4

const MAX_GLOBAL_TRIES = 10_000
const MAX_TRIES_H1 = 200
const MAX_TRIES_H2 = 300


function usage()
    println("Usage: julia --project=. scripts/solve_laman_scip.jl [options]")
    println()
    println("Options:")
    println("  --parity        Run parity presolve before SCIP")
    println("  --no-parity     Solve the generated model directly with SCIP")
    println("  --n N           Number of Laman graph vertices (default: $N)")
    println("  --seed SEED     Random seed (default: $SEED)")
    println("  --R R           Coordinate radius (default: $R)")
    println("  --pH2 P         Henneberg II probability (default: $PH2)")
    println("  --time-limit S  SCIP time limit in seconds, or none (default: unlimited)")
    println("  -h, --help      Show this help")
    return nothing
end


function parse_value(::Type{T}, args::Vector{String}, idx::Int, name::String) where {T}
    idx < length(args) || error("Missing value for $name")
    try
        return parse(T, args[idx + 1])
    catch
        error("Invalid value for $name: $(args[idx + 1])")
    end
end


function parse_optional_float(args::Vector{String}, idx::Int, name::String)
    idx < length(args) || error("Missing value for $name")
    value = args[idx + 1]
    lowercase(value) == "none" && return nothing

    try
        parsed = parse(Float64, value)
        parsed > 0.0 || error("Time limit must be positive, got $value")
        return parsed
    catch err
        err isa ErrorException && rethrow()
        error("Invalid value for $name: $value")
    end
end


function parse_args(args::Vector{String})
    config = (
        n = N,
        radius = R,
        pH2 = PH2,
        seed = SEED,
        use_parity = USE_PARITY_PRESOLVE,
        time_limit_sec = TIME_LIMIT_SEC,
    )

    idx = 1
    while idx <= length(args)
        arg = args[idx]

        if arg == "-h" || arg == "--help"
            usage()
            exit(0)
        elseif arg == "--parity"
            config = merge(config, (; use_parity = true))
        elseif arg == "--no-parity"
            config = merge(config, (; use_parity = false))
        elseif arg == "--n"
            config = merge(config, (; n = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--seed"
            config = merge(config, (; seed = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--R"
            config = merge(config, (; radius = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--pH2"
            config = merge(config, (; pH2 = parse_value(Float64, args, idx, arg)))
            idx += 1
        elseif arg == "--time-limit" || arg == "--time-limit-sec"
            config = merge(config, (; time_limit_sec = parse_optional_float(args, idx, arg)))
            idx += 1
        else
            error("Unknown argument: $arg")
        end

        idx += 1
    end

    return config
end


function build_laman_qp_model(config)
    graph, coords = GE.random_laman_graph(
        config.n;
        R = config.radius,
        pH2 = config.pH2,
        seed = config.seed,
        max_global_tries = MAX_GLOBAL_TRIES,
        max_tries_H1 = MAX_TRIES_H1,
        max_tries_H2 = MAX_TRIES_H2,
    )
    embedded_graph = GE.to_embedded(graph, coords)
    jump_model, _, _ = GE.build_embedding_model(embedded_graph)

    qp_model_builder = QIP.from_moi(backend(jump_model))
    qp_model = QIP.build_model(qp_model_builder)
    PC.fix_vars!(qp_model)

    return graph, qp_model
end


function model_stats(model::PC.QPModel)
    return (
        nvars = length(model.vars),
        nconstraints = length(model.cons),
        infeasible = model.infeasible,
    )
end


function print_model_stats(label::AbstractString, stats)
    println(
        "$label: vars=$(stats.nvars) constraints=$(stats.nconstraints) infeasible=$(stats.infeasible)",
    )
    return nothing
end


function maybe_run_parity_presolve!(model::PC.QPModel, use_parity::Bool)
    use_parity || return nothing

    postsolver = PC.ParityPostsolver(keys(model.vars))
    start_time = time_ns()
    PC.parity_presolve!(model, postsolver)
    elapsed_sec = (time_ns() - start_time) / 1.0e9

    println("parity_presolve_time_sec=$(elapsed_sec)")
    return postsolver
end


function build_scip_model(model::PC.QPModel, time_limit_sec::Union{Nothing, Float64})
    moi_model = QIP.build_moi_model(model, :scip)
    optimizer = getfield(moi_model, :model)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("display/verblevel"), SCIP_VERBLEVEL)
    if time_limit_sec !== nothing
        MOI.set(optimizer, MOI.RawOptimizerAttribute("limits/time"), time_limit_sec)
    end
    return moi_model
end


function variable_id_from_name(name::String)
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
        solution[var_id] = MOI.get(moi_model, MOI.VariablePrimal(), variable)
    end

    return solution
end


function print_solution_sample(label::AbstractString, solution::Dict{PC.VarId, Float64}; limit::Int = 12)
    println("$label: $(length(solution)) values")
    isempty(solution) && return nothing

    shown = 0
    for var_id in sort!(collect(keys(solution)))
        shown += 1
        print("x$var_id=$(solution[var_id])")
        if shown == min(limit, length(solution))
            println()
            break
        end
        print(", ")
    end

    length(solution) > limit && println("... $(length(solution) - limit) more values omitted")
    return nothing
end


function has_feasible_solution(primal_status)
    return primal_status == MOI.FEASIBLE_POINT || primal_status == MOI.NEARLY_FEASIBLE_POINT
end


function solve_with_scip(model::PC.QPModel, time_limit_sec::Union{Nothing, Float64})
    moi_model = build_scip_model(model, time_limit_sec)

    start_time = time_ns()
    MOI.optimize!(moi_model)
    solve_time_sec = (time_ns() - start_time) / 1.0e9

    result_count = MOI.get(moi_model, MOI.ResultCount())
    termination_status = MOI.get(moi_model, MOI.TerminationStatus())
    primal_status = result_count > 0 ? MOI.get(moi_model, MOI.PrimalStatus()) : MOI.NO_SOLUTION

    return (
        moi_model = moi_model,
        solve_time_sec = solve_time_sec,
        result_count = result_count,
        termination_status = termination_status,
        primal_status = primal_status,
    )
end


function main(args::Vector{String})
    config = parse_args(args)

    graph, qp_model = build_laman_qp_model(config)
    original_stats = model_stats(qp_model)

    println("n=$(nv(graph)) edges=$(ne(graph)) seed=$(config.seed) R=$(config.radius) pH2=$(config.pH2)")
    println("use_parity_presolve=$(config.use_parity)")
    println("scip_time_limit_sec=$(config.time_limit_sec === nothing ? "unlimited" : string(config.time_limit_sec))")
    println("scip_display_verblevel=$SCIP_VERBLEVEL")
    print_model_stats("before_parity", original_stats)

    postsolver = maybe_run_parity_presolve!(qp_model, config.use_parity)
    presolved_stats = model_stats(qp_model)
    print_model_stats("after_parity", presolved_stats)

    if qp_model.infeasible
        println("SCIP skipped because parity presolve proved the model infeasible.")
        return nothing
    end

    result = solve_with_scip(qp_model, config.time_limit_sec)
    println("scip_termination_status=$(result.termination_status)")
    println("scip_primal_status=$(result.primal_status)")
    println("scip_result_count=$(result.result_count)")
    println("scip_solve_time_sec=$(result.solve_time_sec)")

    if result.result_count > 0 && has_feasible_solution(result.primal_status)
        reduced = reduced_solution(result.moi_model)
        print_solution_sample("reduced_solution_sample", reduced)

        if postsolver !== nothing
            original = PC.postsolve(postsolver, reduced)
            print_solution_sample("original_solution_sample", original)
        end
    end

    return nothing
end


main(copy(ARGS))
