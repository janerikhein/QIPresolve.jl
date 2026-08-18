#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuMP: backend
using Printf

import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC
using QIPresolve.InstanceGeneration: generate_random_qip_model

const DEFAULT_N_VARS = 100
const DEFAULT_N_CONSTRAINTS = 200
const DEFAULT_SEED = 42

const RANDOM_QIP_KWARGS = (
    p_con_eq = 0.0,
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
    constraint_slack_range = -10:10,
)

log_domain_sum(model::PC.QPModel) = sum(log(var.ub - var.lb + 1) for var in values(model.vars); init = 0.0)
domain_sum(model::PC.QPModel) = sum(var.ub - var.lb + 1 for var in values(model.vars); init = 0.0)

function avg_constraint_bound_range(model::PC.QPModel)
    isempty(model.cons) && return 0.0
    return sum(con.rhs - con.lhs for con in model.cons; init = 0.0) / length(model.cons)
end

function parse_int_arg(args::Vector{String}, idx::Int, name::String, default::Int)
    length(args) < idx && return default
    try
        return parse(Int, args[idx])
    catch
        error("Invalid $name: $(args[idx])")
    end
end

function usage()
    println("Usage: julia --project=. scripts/presolve_random_qip_metrics.jl [nvars] [ncons] [seed]")
    return println("Defaults: nvars=$DEFAULT_N_VARS ncons=$DEFAULT_N_CONSTRAINTS seed=$DEFAULT_SEED")
end

function build_random_qip_core_model(nvars::Int, ncons::Int, seed::Int)
    jump_model, _ = generate_random_qip_model(nvars, ncons; RANDOM_QIP_KWARGS..., seed = seed)
    return QIP.build_model(QIP.from_moi(backend(jump_model)))
end

function metrics(model::PC.QPModel)
    return (
        log_domain_sum = log_domain_sum(model),
        domain_sum = domain_sum(model),
        avg_constraint_bound_range = avg_constraint_bound_range(model),
    )
end

function print_config(nvars::Int, ncons::Int, seed::Int)
    println("Random QIP presolve metrics")
    println("nvars = $nvars")
    println("nconstraints = $ncons")
    println("seed = $seed")
    println("residue_strategy = $(QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY)")
    println("residue_threshold = $(QIP.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD)")
    println("treewidth_threshold = $(QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD)")
    println("p_var_bilin = $(RANDOM_QIP_KWARGS.p_var_bilin)")
    println()
    return nothing
end

function print_metrics_table(before, after)
    @printf("%-32s %14s %14s\n", "metric", "before", "after")
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

    length(args) <= 3 || error("Expected at most 3 positional arguments, got $(length(args))")

    nvars = parse_int_arg(args, 1, "nvars", DEFAULT_N_VARS)
    ncons = parse_int_arg(args, 2, "ncons", DEFAULT_N_CONSTRAINTS)
    seed = parse_int_arg(args, 3, "seed", DEFAULT_SEED)

    model = build_random_qip_core_model(nvars, ncons, seed)
    infeasible_before = model.infeasible
    before = metrics(model)

    result = QIP.presolve!(model)
    presolved_model = result.model
    after = metrics(presolved_model)

    print_config(nvars, ncons, seed)
    println("infeasible_before = $infeasible_before")
    println("infeasible_after = $(presolved_model.infeasible)")
    println()
    print_metrics_table(before, after)

    return nothing
end

main(ARGS)
