using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using JuMP: backend
using Printf

import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF
import QIPresolve as QIP
import SCIP

PC = QIP.PresolvingCore
GE = QIP.InstanceGeneration

const N_MIN = 10
const N_MAX = 150
const N_STEP = 10
const N_INSTANCES = 1
const R = 20
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


function write_stats_row(io, nvertices::Int, instance_idx::Int, seed::Int, stage::String, stats, error)
    error_text = error === nothing ? "" : sprint(showerror, error)
    @printf(
        io,
        "%d,%d,%d,%s,%d,%d,%d,%d,%.12g,%.12g,%s\n",
        nvertices,
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


function write_error_rows(io, nvertices::Int, instance_idx::Int, seed::Int, stages, error)
    empty_stats = (
        nvars = 0,
        total_constraints = 0,
        quadratic_constraints = 0,
        linear_constraints = 0,
        log_domain_sum = NaN,
        total_domain_sum = NaN,
    )

    for stage in stages
        write_stats_row(io, nvertices, instance_idx, seed, stage, empty_stats, error)
    end

    return nothing
end


function run_instance(io, nvertices::Int, instance_idx::Int)
    seed = SEED_BASE + 10_000 * nvertices + instance_idx
    stages = ("original", "after_scip_1", "after_parity", "after_scip_2")

    qp_model = try
        build_laman_qp_model(seed, nvertices)
    catch err
        write_error_rows(io, nvertices, instance_idx, seed, stages, err)
        return false
    end

    write_stats_row(io, nvertices, instance_idx, seed, "original", model_stats(qp_model), nothing)

    mktempdir() do temp_dir
        scip_model = try
            run_scip_presolve(qp_model, temp_dir, "after_scip_1")
        catch err
            write_error_rows(io, nvertices, instance_idx, seed, stages[2:end], err)
            return false
        end
        write_stats_row(io, nvertices, instance_idx, seed, "after_scip_1", model_stats(scip_model), nothing)

        parity_model = try
            model_copy = deepcopy(scip_model)
            run_parity_presolve!(model_copy)
        catch err
            write_error_rows(io, nvertices, instance_idx, seed, stages[3:end], err)
            return false
        end
        write_stats_row(io, nvertices, instance_idx, seed, "after_parity", model_stats(parity_model), nothing)

        scip_model_2 = try
            run_scip_presolve(parity_model, temp_dir, "after_scip_2")
        catch err
            write_error_rows(io, nvertices, instance_idx, seed, stages[4:end], err)
            return false
        end
        write_stats_row(io, nvertices, instance_idx, seed, "after_scip_2", model_stats(scip_model_2), nothing)
    end

    return true
end


function run_sweep(nmin::Int, nmax::Int, step::Int, ninstances::Int)
    nvalues = sweep_nvalues(nmin, nmax, step)
    results_dir = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(results_dir)

    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_path = joinpath(results_dir, "$(splitext(basename(@__FILE__))[1])-$timestamp.csv")
    config_line =
        "sweep_config: n_min=$nmin n_max=$nmax step=$step ninstances=$ninstances R=$R pH2=$PH2 seed_base=$SEED_BASE"
    header =
        "nvertices,instance_idx,seed,stage,nvars,total_constraints,quadratic_constraints,linear_constraints,log_domain_sum,total_domain_sum,error"

    println(config_line)
    println(header)

    open(output_path, "w") do io
        println(io, "# ", config_line)
        println(io, header)

        for nvertices in nvalues
            for instance_idx in 1:ninstances
                ok = run_instance(io, nvertices, instance_idx)
                status = ok ? "ok" : "error"
                println("completed nvertices=$nvertices instance_idx=$instance_idx status=$status")
                flush(io)
            end
        end
    end

    println("wrote results to $output_path")
    return output_path
end


run_sweep(N_MIN, N_MAX, N_STEP, N_INSTANCES)
