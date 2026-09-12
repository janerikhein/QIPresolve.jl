#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

module PresolveLpScipStatsScript

import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
import QIPresolve as QIP
import SCIP

const SCIP_TIME_LIMIT_SEC = 3600.0
const USAGE = "Usage: julia --project=. scripts/presolve_lp_scip_stats.jl instance.lp"

seconds_since(start_time) = (time_ns() - start_time) / 1.0e9

function read_lp(path::AbstractString)
    isfile(path) || error("LP file not found: $path")
    model = FF.Model(format = FF.FORMAT_LP)
    MOI.read_from_file(model, path)

    integer_variables = Set{MOI.VariableIndex}()
    for set_type in (MOI.Integer, MOI.ZeroOne)
        for ci in MOI.get(model, MOI.ListOfConstraintIndices{MOI.VariableIndex, set_type}())
            push!(integer_variables, MOI.get(model, MOI.ConstraintFunction(), ci))
        end
    end
    all(in(integer_variables), MOI.get(model, MOI.ListOfVariableIndices())) ||
        error("Unsupported LP: QIPresolve supports only integer and binary variables")
    MOI.get(model, MOI.ObjectiveFunctionType()) in
        (MOI.VariableIndex, MOI.ScalarAffineFunction{Float64}) ||
        error("Unsupported LP: this script supports affine objectives only")
    return model
end

function build_core_model(lp_model)
    builder = QIP.from_moi(lp_model)
    # The core stores domains as bounds, so retain the bounds implied by Binary.
    for (id, info) in builder.vars
        info.var_type == :bin || continue
        QIP.ModelIO.register_var_info!(builder, id; lb = 0.0, ub = 1.0)
    end
    return QIP.build_model(builder)
end

function configure_scip!(optimizer::SCIP.Optimizer)
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("timing/clocktype"), 2) # wall clock
    MOI.set(optimizer, MOI.TimeLimitSec(), SCIP_TIME_LIMIT_SEC)
    return optimizer
end

function scip_metrics(optimizer::SCIP.Optimizer)
    scip = optimizer.inner
    stage = SCIP.SCIPgetStage(scip)
    # A time limit can stop SCIP before its presolving clock is available.
    presolve_time = if SCIP.SCIP_STAGE_INITPRESOLVE <= stage <= SCIP.SCIP_STAGE_SOLVED
        SCIP.SCIPgetPresolvingTime(scip)
    else
        0.0
    end
    # SCIP's solving clock includes its presolving time.
    return (
        scip_presolve_time_sec = presolve_time,
        scip_solving_time_sec = max(0.0, SCIP.SCIPgetSolvingTime(scip) - presolve_time),
        status = MOI.get(optimizer, MOI.RawStatusString()),
    )
end

function solve_with_scip(source, start_time, qip_presolve_time_sec)
    optimizer = SCIP.Optimizer()
    try
        configure_scip!(optimizer)
        model = MOI.Bridges.full_bridge_optimizer(optimizer, Float64)
        MOI.copy_to(model, source)
        MOI.optimize!(model)
        wall_time_sec = seconds_since(start_time)
        return (; wall_time_sec, qip_presolve_time_sec, scip_metrics(optimizer)...)
    finally
        SCIP.free_scip(optimizer.inner)
    end
end

function run_experiment(lp_model; presolve::Bool)
    start_time = time_ns()
    qip_presolve_time_sec = 0.0
    source = lp_model
    if presolve
        model = build_core_model(lp_model)
        presolve_start = time_ns()
        QIP.presolve!(model; enable_parity = true, enable_residue = true)
        qip_presolve_time_sec = seconds_since(presolve_start)
        if model.infeasible
            return (
                wall_time_sec = seconds_since(start_time),
                qip_presolve_time_sec = qip_presolve_time_sec,
                scip_presolve_time_sec = 0.0,
                scip_solving_time_sec = 0.0,
                status = "INFEASIBLE (QIPresolve; SCIP skipped)",
            )
        end
        source = QIP.build_moi_model(model)
    end
    return solve_with_scip(source, start_time, qip_presolve_time_sec)
end

function print_result(io::IO, label::AbstractString, result)
    println(io, label)
    println(io, "  Full wall time (s): ", result.wall_time_sec)
    println(io, "  Parity/residue presolve time (s): ", result.qip_presolve_time_sec)
    println(io, "  SCIP presolve time (s): ", result.scip_presolve_time_sec)
    println(io, "  SCIP solving time (s): ", result.scip_solving_time_sec)
    println(io, "  SCIP status: ", result.status)
end

"""
Run both experiments and print their timings in seconds and SCIP statuses.
Wall time includes model preparation, reductions, and optimization, excluding
shared LP loading, Julia startup, reporting, and solver cleanup. Julia compilation
within a measured call is included; no warm-up runs are performed.
"""
function main(args::Vector{String} = copy(ARGS); io::IO = stdout)
    if length(args) == 1 && only(args) in ("-h", "--help")
        println(io, USAGE)
        println(io, "Compare SCIP alone with parity/residue presolve followed by SCIP.")
        println(io, "Each SCIP run has a 1800-second wall-clock limit. Results are printed in seconds.")
        return nothing
    end
    length(args) == 1 && !startswith(only(args), "-") || error(USAGE)
    lp_model = read_lp(only(args))
    original = run_experiment(lp_model; presolve = false)
    print_result(io, "SCIP only", original)
    println(io)
    presolved = run_experiment(lp_model; presolve = true)
    print_result(io, "Parity/residue + SCIP", presolved)
    return (; original, presolved)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    PresolveLpScipStatsScript.main()
end
