using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Printf
using Random

import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
import QIPresolve as QIP
import SCIP

PC = QIP.PresolvingCore

const QuadTerm = Tuple{Float64, Int, Int}
const LinTerm = Tuple{Float64, Int}

const N_VARS = 200
const N_INSTANCES = 1
const CONSTRAINT_STEP_FACTOR = 0.1
const MAX_CONSTRAINT_FACTOR = 10
const VAR_LB = 0
const VAR_UB = 100
const SEED_BASE = 19
const MAX_MODEL_TRIES = 100
const DENSITY = 0.05


function random_quadratic_terms(
        rng::AbstractRNG, nvars::Int; force_odd_diagonal::Bool = false
    )::Vector{QuadTerm}
    terms = QuadTerm[]

    for i in 1:nvars
        rand(rng) < DENSITY || continue
        coeff = rand(rng, -500:500)
        coeff == 0 && continue
        push!(terms, (Float64(coeff), i, i))
    end

    for i in 1:(nvars - 1)
        for j in (i + 1):nvars
            rand(rng) < DENSITY || continue
            coeff = rand(rng, -500:500) * 2
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


function build_random_qp_model(seed::Int, nvars::Int, ncons::Int)
    rng = MersenneTwister(seed)
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
        parity_model = PC.build_parity_model(model)
        if !isempty(parity_model.pos_to_var_id)
            return model
        end
    end

    error(
        "Failed to generate a model with nontrivial parity structure after $MAX_MODEL_TRIES attempts."
    )
end


function has_nonzero_quadratic_terms(expr::PC.QuadExpr)
    var_ids = sort!(collect(PC.vars(expr)))
    for (i, var_id_1) in enumerate(var_ids)
        for var_id_2 in @view var_ids[i:end]
            isapprox(PC.get_quad_coeff(expr, var_id_1, var_id_2), 0.0) || return true
        end
    end

    return false
end


function model_stats(model::PC.QPModel)
    total_constraints = length(model.cons)
    quadratic_constraints = count(con -> has_nonzero_quadratic_terms(con.qe), model.cons)
    linear_constraints = total_constraints - quadratic_constraints
    log_total = sum(log(var.ub - var.lb + 1) for var in values(model.vars); init = 0.0)
    domain_total = sum(var.ub - var.lb + 1 for var in values(model.vars); init = 0.0)

    return (
        nvars = length(model.vars),
        total_constraints = total_constraints,
        quadratic_constraints = quadratic_constraints,
        linear_constraints = linear_constraints,
        log_domain_sum = log_total,
        total_domain_sum = domain_total,
    )
end


function run_parity_presolve!(model::PC.QPModel)
    propagator = PC.PropagationManager(Int[])

    while !model.infeasible
        stats = PC.parity_presolve_phase!(model, propagator)
        PC.scale_constraints_gcd!(model)
        stats.changed || break
    end

    return model
end


function coerce_integral_continuous_variables!(builder)
    for (var_id, info) in builder.vars
        info.var_type == :cont || continue
        isfinite(info.lb) && isfinite(info.ub) || continue
        isinteger(info.lb) && isinteger(info.ub) || continue

        var_type = info.lb == 0.0 && info.ub == 1.0 ? :bin : :int
        builder.vars[var_id] = QIP.ModelIO.VarInfo(info.lb, info.ub, var_type)
    end

    return builder
end


function read_lp_as_qp_model(path::AbstractString)
    moi_model = FF.Model(format = FF.FORMAT_LP)
    MOI.read_from_file(moi_model, path)

    qp_model_builder = QIP.from_moi(moi_model)
    coerce_integral_continuous_variables!(qp_model_builder)
    return QIP.build_model(qp_model_builder)
end


function run_scip_presolve(qp_model::PC.QPModel, temp_dir::AbstractString, label::AbstractString)
    moi_model = QIP.build_moi_model(qp_model, :scip)
    scip_backend = getfield(moi_model, :model)
    MOI.set(scip_backend, MOI.RawOptimizerAttribute("display/verblevel"), 0)

    scip = scip_backend.inner
    lp_path = joinpath(temp_dir, "$label.lp")
    SCIP.@SCIP_CALL SCIP.SCIPpresolve(scip)
    SCIP.@SCIP_CALL SCIP.SCIPwriteTransProblem(scip, lp_path, "lp", true)

    return read_lp_as_qp_model(lp_path)
end


function csv_escape(value)
    text = replace(string(value), '\n' => ' ', '\r' => ' ')
    if occursin(',', text) || occursin('"', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end

    return text
end


function write_stats_row(io, nvars::Int, ncons::Int, instance_idx::Int, seed::Int, stage::String, stats, error)
    error_text = error === nothing ? "" : sprint(showerror, error)
    @printf(
        io,
        "%d,%d,%d,%d,%s,%d,%d,%d,%d,%.12g,%.12g,%s\n",
        nvars,
        ncons,
        instance_idx,
        seed,
        stage,
        stats.nvars,
        stats.total_constraints,
        stats.quadratic_constraints,
        stats.linear_constraints,
        stats.log_domain_sum,
        stats.total_domain_sum,
        csv_escape(error_text),
    )
    return nothing
end


function write_error_rows(io, nvars::Int, ncons::Int, instance_idx::Int, seed::Int, stages, error)
    empty_stats = (
        nvars = 0,
        total_constraints = 0,
        quadratic_constraints = 0,
        linear_constraints = 0,
        log_domain_sum = NaN,
        total_domain_sum = NaN,
    )

    for stage in stages
        write_stats_row(io, nvars, ncons, instance_idx, seed, stage, empty_stats, error)
    end

    return nothing
end


function run_instance(io, nvars::Int, ncons::Int, instance_idx::Int)
    seed = SEED_BASE + 10_000 * ncons + instance_idx
    stages = ("original", "after_scip_1", "after_parity", "after_scip_2")

    qp_model = try
        build_random_qp_model(seed, nvars, ncons)
    catch err
        write_error_rows(io, nvars, ncons, instance_idx, seed, stages, err)
        return false
    end

    write_stats_row(io, nvars, ncons, instance_idx, seed, "original", model_stats(qp_model), nothing)

    mktempdir() do temp_dir
        scip_model = try
            run_scip_presolve(qp_model, temp_dir, "after_scip_1")
        catch err
            write_error_rows(io, nvars, ncons, instance_idx, seed, stages[2:end], err)
            return false
        end
        write_stats_row(io, nvars, ncons, instance_idx, seed, "after_scip_1", model_stats(scip_model), nothing)

        parity_model = try
            model_copy = deepcopy(scip_model)
            run_parity_presolve!(model_copy)
        catch err
            write_error_rows(io, nvars, ncons, instance_idx, seed, stages[3:end], err)
            return false
        end
        write_stats_row(io, nvars, ncons, instance_idx, seed, "after_parity", model_stats(parity_model), nothing)

        scip_model_2 = try
            run_scip_presolve(parity_model, temp_dir, "after_scip_2")
        catch err
            write_error_rows(io, nvars, ncons, instance_idx, seed, stages[4:end], err)
            return false
        end
        write_stats_row(io, nvars, ncons, instance_idx, seed, "after_scip_2", model_stats(scip_model_2), nothing)
    end

    return true
end


function run_sweep(nvars::Int, ninstances::Int)
    step = max(1, round(Int, CONSTRAINT_STEP_FACTOR * nvars))
    max_cons = round(Int, MAX_CONSTRAINT_FACTOR * nvars)
    results_dir = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(results_dir)

    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_path = joinpath(results_dir, "$(splitext(basename(@__FILE__))[1])-$timestamp.csv")
    config_line = "sweep_config: nvars=$nvars ninstances=$ninstances step=$step max_constraints=$max_cons seed_base=$SEED_BASE density=$DENSITY"
    header =
        "nvars,nconstraints,instance_idx,seed,stage,remaining_vars,total_constraints,quadratic_constraints,linear_constraints,log_domain_sum,total_domain_sum,error"

    println(config_line)
    println(header)

    open(output_path, "w") do io
        println(io, "# ", config_line)
        println(io, header)

        for ncons in step:step:max_cons
            for instance_idx in 1:ninstances
                ok = run_instance(io, nvars, ncons, instance_idx)
                status = ok ? "ok" : "error"
                println("completed nconstraints=$ncons instance_idx=$instance_idx status=$status")
                flush(io)
            end
        end
    end

    println("wrote results to $output_path")
    return output_path
end


run_sweep(N_VARS, N_INSTANCES)
