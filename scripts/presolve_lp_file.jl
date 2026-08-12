#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
using Printf

import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC

log_domain_sum(model::PC.QPModel) = sum(log(var.ub - var.lb + 1) for var in values(model.vars); init = 0.0)
domain_sum(model::PC.QPModel) = sum(var.ub - var.lb + 1 for var in values(model.vars); init = 0.0)

function avg_constraint_bound_range(model::PC.QPModel)
    isempty(model.cons) && return 0.0
    return sum(con.rhs - con.lhs for con in model.cons; init = 0.0) / length(model.cons)
end

function usage()
    return println("Usage: julia --project=. scripts/presolve_lp_file.jl path/to/model.lp")
end

function load_lp_core_model(file_path::AbstractString)
    moi_model = FF.Model(format = FF.FORMAT_LP)
    MOI.read_from_file(moi_model, file_path)
    return QIP.build_model(QIP.from_moi(moi_model))
end

function metrics(model::PC.QPModel)
    return (
        nvars = length(model.vars),
        nconstraints = length(model.cons),
        log_domain_sum = log_domain_sum(model),
        domain_sum = domain_sum(model),
        avg_constraint_bound_range = avg_constraint_bound_range(model),
    )
end

function print_config(file_path::AbstractString)
    println("LP file presolve metrics")
    println("input_path = $file_path")
    println("residue_strategy = $(QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY)")
    println("residue_threshold = $(QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD)")
    println("treewidth_threshold = $(QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD)")
    println()
    return nothing
end

function print_metrics_table(before, after)
    @printf("%-32s %14s %14s\n", "metric", "before", "after")
    @printf("%-32s %14d %14d\n", "nvars", before.nvars, after.nvars)
    @printf("%-32s %14d %14d\n", "nconstraints", before.nconstraints, after.nconstraints)
    @printf("%-32s %14.6f %14.6f\n", "log_domain_sum", before.log_domain_sum, after.log_domain_sum)
    @printf("%-32s %14.6f %14.6f\n", "domain_sum", before.domain_sum, after.domain_sum)
    @printf(
        "%-32s %14.6f %14.6f\n",
        "avg_constraint_bound_range",
        before.avg_constraint_bound_range,
        after.avg_constraint_bound_range,
    )
    return nothing
end

function main(args::Vector{String})
    if !isempty(args) && args[1] in ("-h", "--help")
        usage()
        return nothing
    end

    length(args) == 1 || error("Expected exactly one LP file path, got $(length(args))")

    file_path = args[1]
    isfile(file_path) || error("LP file not found: $file_path")

    model = load_lp_core_model(file_path)
    infeasible_before = model.infeasible
    before = metrics(model)

    result = QIP.presolve!(model)
    presolved_model = result.model
    after = metrics(presolved_model)

    print_config(file_path)
    println("infeasible_before = $infeasible_before")
    println("infeasible_after = $(presolved_model.infeasible)")
    println()
    print_metrics_table(before, after)

    return nothing
end

main(ARGS)
