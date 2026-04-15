using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Printf
using Random

import QIPresolve as QIP

PC = QIP.PresolvingCore

const QuadTerm = Tuple{Float64, Int, Int}
const LinTerm = Tuple{Float64, Int}

const N_VARS = 100
const N_INSTANCES = 50
const CONSTRAINT_STEP_FACTOR = 0.1
const MAX_CONSTRAINT_FACTOR = 5
const VAR_LB = 0
const VAR_UB = 100
const SEED_BASE = 19
const MAX_MODEL_TRIES = 100
const DENSITY = 0.005


function random_quadratic_terms(
        rng::AbstractRNG, nvars::Int; force_odd_diagonal::Bool = false
    )::Vector{QuadTerm}
    terms = QuadTerm[]

    for i in 1:nvars
        rand(rng) < DENSITY || continue
        coeff = rand(rng, -5:5)
        coeff == 0 && continue
        push!(terms, (Float64(coeff), i, i))
    end

    for i in 1:(nvars - 1)
        for j in (i + 1):nvars
            rand(rng) < DENSITY || continue
            coeff = rand(rng, -5:5) * 2
            coeff == 0 && continue
            push!(terms, (Float64(coeff), i, j))
        end
    end

    if force_odd_diagonal
        has_odd_diag = any(term -> term[2] == term[3] && isodd(round(Int, term[1])), terms)
        if !has_odd_diag
            vid = rand(rng, 1:nvars)
            replaced = false
            for idx in eachindex(terms)
                coeff, i, j = terms[idx]
                if i == vid && j == vid
                    terms[idx] = (coeff >= 0 ? 1.0 : -1.0, vid, vid)
                    replaced = true
                    break
                end
            end
            replaced || push!(terms, (rand(rng, Bool) ? 1.0 : -1.0, vid, vid))
        end
    end

    isempty(terms) && return random_quadratic_terms(
        rng, nvars; force_odd_diagonal = force_odd_diagonal
    )
    return terms
end


function build_random_constraint(
        rng::AbstractRNG, con_id::Int, x_star::Vector{Int}; force_odd_diagonal::Bool = false
    )::PC.Constraint
    quad_terms = random_quadratic_terms(rng, length(x_star); force_odd_diagonal = force_odd_diagonal)
    qe = PC.QuadExpr(quad_terms, LinTerm[]; constant = 0.0)
    rhs = Float64(PC.eval_full(qe, x_star))
    return PC.Constraint(con_id, qe, rhs, rhs)
end


function build_random_qp_model(rng::AbstractRNG, nvars::Int, ncons::Int)
    x_star = rand(rng, VAR_LB:VAR_UB, nvars)

    vars = Dict(i => PC.IntVar(Float64(VAR_LB), Float64(VAR_UB)) for i in 1:nvars)
    obj = PC.QuadExpr(QuadTerm[], LinTerm[]; constant = 0.0)

    for _ in 1:MAX_MODEL_TRIES
        cons = Vector{PC.Constraint}(undef, ncons)
        for con_id in 1:ncons
            cons[con_id] = build_random_constraint(
                rng, con_id, x_star; force_odd_diagonal = con_id == 1
            )
        end

        model = PC.QPModel(copy(vars), cons, obj, :min)
        xor_model = PC.build_parity_model(model)
        if !isempty(xor_model.pos_to_var_id)
            return model
        end
    end

    error(
        "Failed to generate a model with nontrivial parity structure after $MAX_MODEL_TRIES attempts."
    )
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


function evaluate_instance(rng, nvars::Int, ncons::Int)
    build_error = nothing
    qp_model = nothing

    try
        qp_model = build_random_qp_model(rng, nvars, ncons)
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


function run_sweep(nvars::Int, ninstances::Int)
    step = max(1, round(Int, CONSTRAINT_STEP_FACTOR * nvars))
    max_cons = round(Int, MAX_CONSTRAINT_FACTOR * nvars)
    config_line = "sweep_config: nvars=$nvars ninstances=$ninstances step=$step max_constraints=$max_cons seed_base=$SEED_BASE"
    console_header =
        "m\tavg_domain_before\tavg_domain_after\tavg_log_domain_before\tavg_log_domain_after\tavg_presolve_runtime_sec\tsolved_fraction\tn_instances\tn_errors"
    csv_header =
        "m,avg_domain_before,avg_domain_after,avg_log_domain_before,avg_log_domain_after,avg_presolve_runtime_sec,solved_fraction,n_instances,n_errors"
    results_dir = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(results_dir)
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_path = joinpath(results_dir, "$(splitext(basename(@__FILE__))[1])-$timestamp.csv")

    println(config_line)
    println(console_header)

    open(output_path, "w") do output_io
        println(output_io, "# ", config_line)
        println(output_io, csv_header)

        termination_reason = "stopped because m reached 10n"
        rng = MersenneTwister(SEED_BASE)

        for ncons in step:step:max_cons
            total_domain_before = 0.0
            total_domain_after = 0.0
            total_before = 0.0
            total_after = 0.0
            total_presolve_runtime = 0.0
            solved_count = 0
            error_count = 0
            successful_count = 0

            for instance_idx in 1:ninstances
                result = evaluate_instance(rng, nvars, ncons)
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

            @printf(
                "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.3f\t%d\t%d\n",
                ncons,
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
                "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.3f,%d,%d\n",
                ncons,
                avg_domain_before,
                avg_domain_after,
                avg_before,
                avg_after,
                avg_presolve_runtime,
                solved_fraction,
                ninstances,
                error_count,
            )

            if solved_fraction == 1.0
                termination_reason = "stopped because solved_fraction == 1.0"
                break
            end
        end

        println(termination_reason)
    end
    return nothing
end


run_sweep(N_VARS, N_INSTANCES)
