using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Printf
using JuMP: backend

import QIPresolve as QIP

PC = QIP.PresolvingCore
GE = QIP.GraphEmbedding

const N_MIN = 3
const N_MAX = 150
const N_STEP = 5
const N_INSTANCES = 50
const R = 100
const PH2 = 0.5
const SEED_BASE = 17
const MAX_GLOBAL_TRIES = 10_000
const MAX_TRIES_H1 = 200
const MAX_TRIES_H2 = 300


function sweep_nvalues(nmin::Int, nmax::Int, step::Int)
    nvalues = collect(nmin:step:nmax)
    isempty(nvalues) && return [nmax]
    last(nvalues) == nmax && return nvalues
    push!(nvalues, nmax)
    return nvalues
end


function build_laman_qp_model(seed::Int, nvertices::Int)
    graph, coords = GE.random_laman_graph(
        nvertices;
        R = R,
        pH2 = PH2,
        seed = seed,
        max_global_tries = MAX_GLOBAL_TRIES,
        max_tries_H1 = MAX_TRIES_H1,
        max_tries_H2 = MAX_TRIES_H2,
    )
    emb_graph = GE.to_embedded(graph, coords)
    jump_model, _, _ = GE.build_embedding_model(emb_graph)

    qp_model_builder = QIP.from_moi(backend(jump_model))
    qp_model = QIP.build_model(qp_model_builder)
    PC.fix_vars!(qp_model)
    return qp_model
end


log_domain_sum(model::PC.QPModel) = sum(log(var.ub - var.lb + 1) for var in values(model.vars); init = 0.0)
domain_sum(model::PC.QPModel) = sum(var.ub - var.lb + 1 for var in values(model.vars); init = 0.0)


function run_presolve!(qp_model::PC.QPModel)
    propagator = PC.PropagationManager(Int[])

    while !qp_model.infeasible
        stats = PC.parity_presolve_phase!(qp_model, propagator)
        PC.scale_constraints_gcd!(qp_model)
        stats.changed || break
    end

    return nothing
end


function evaluate_instance(seed::Int, nvertices::Int)
    build_error = nothing
    qp_model = nothing

    try
        qp_model = build_laman_qp_model(seed, nvertices)
    catch err
        build_error = err
    end

    build_error === nothing || return (
        domain_before = NaN,
        domain_after = NaN,
        log_before = NaN,
        log_after = NaN,
        presolve_runtime_sec = NaN,
        solved = false,
        error = build_error,
    )

    domain_before = domain_sum(qp_model)
    log_before = log_domain_sum(qp_model)
    run_error = nothing
    presolve_runtime_sec = NaN
    start_time_ns = time_ns()

    try
        run_presolve!(qp_model)
        presolve_runtime_sec = (time_ns() - start_time_ns) / 1.0e9
    catch err
        presolve_runtime_sec = (time_ns() - start_time_ns) / 1.0e9
        run_error = err
    end

    domain_after = run_error === nothing ? domain_sum(qp_model) : NaN
    log_after = run_error === nothing ? log_domain_sum(qp_model) : NaN
    solved = run_error === nothing && !qp_model.infeasible && isempty(qp_model.vars)

    return (
        domain_before = domain_before,
        domain_after = domain_after,
        log_before = log_before,
        log_after = log_after,
        presolve_runtime_sec = presolve_runtime_sec,
        solved = solved,
        error = run_error,
    )
end


function run_sweep(nmin::Int, nmax::Int, step::Int, ninstances::Int)
    nvalues = sweep_nvalues(nmin, nmax, step)
    config_line =
        "sweep_config: n_min=$nmin n_max=$nmax step=$step ninstances=$ninstances R=$R pH2=$PH2 seed_base=$SEED_BASE"
    console_header =
        "nvertices\ttotal_edges\tavg_domain_before\tavg_domain_after\tavg_log_domain_before\tavg_log_domain_after\tavg_presolve_runtime_sec\tsolved_fraction\tn_instances\tn_errors"
    csv_header =
        "nvertices,total_edges,avg_domain_before,avg_domain_after,avg_log_domain_before,avg_log_domain_after,avg_presolve_runtime_sec,solved_fraction,n_instances,n_errors"
    results_dir = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(results_dir)
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_path = joinpath(results_dir, "$(splitext(basename(@__FILE__))[1])-$timestamp.csv")

    println(config_line)
    println(console_header)

    open(output_path, "w") do output_io
        println(output_io, "# ", config_line)
        println(output_io, csv_header)

        for nvertices in nvalues
            total_domain_before = 0.0
            total_domain_after = 0.0
            total_before = 0.0
            total_after = 0.0
            total_presolve_runtime = 0.0
            solved_count = 0
            error_count = 0
            successful_count = 0

            for instance_idx in 1:ninstances
                seed = SEED_BASE + 10_000 * nvertices + instance_idx
                result = evaluate_instance(seed, nvertices)

                if result.error === nothing
                    total_domain_before += result.domain_before
                    total_domain_after += result.domain_after
                    total_before += result.log_before
                    total_after += result.log_after
                    total_presolve_runtime += result.presolve_runtime_sec
                    successful_count += 1
                else
                    error_count += 1
                end

                solved_count += result.solved ? 1 : 0
            end

            avg_domain_before = successful_count == 0 ? NaN : total_domain_before / successful_count
            avg_domain_after = successful_count == 0 ? NaN : total_domain_after / successful_count
            avg_before = successful_count == 0 ? NaN : total_before / successful_count
            avg_after = successful_count == 0 ? NaN : total_after / successful_count
            avg_presolve_runtime = successful_count == 0 ? NaN : total_presolve_runtime / successful_count
            solved_fraction = solved_count / ninstances
            total_edges = 2 * nvertices - 3

            @printf(
                "%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.3f\t%d\t%d\n",
                nvertices,
                total_edges,
                avg_domain_before,
                avg_domain_after,
                avg_before,
                avg_after,
                avg_presolve_runtime,
                solved_fraction,
                ninstances,
                error_count,
            )
            @printf(
                output_io,
                "%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.3f,%d,%d\n",
                nvertices,
                total_edges,
                avg_domain_before,
                avg_domain_after,
                avg_before,
                avg_after,
                avg_presolve_runtime,
                solved_fraction,
                ninstances,
                error_count,
            )
        end
    end

    println("completed sweep because nvertices reached n_max=$nmax")
    return nothing
end


run_sweep(N_MIN, N_MAX, N_STEP, N_INSTANCES)
