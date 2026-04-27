#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
import MathOptInterface as MOI
import QIPresolve as QIP

PC = QIP.PresolvingCore

const QuadTerm = Tuple{Float64, Int, Int}
const LinTerm = Tuple{Float64, Int}

const N_VARS = 20
const N_CONSTRAINTS = 20
const VAR_LB = 0
const VAR_UB = 100
const SEED = 19
const DENSITY = 0.05
const USE_PARITY_PRESOLVE = false
const TIME_LIMIT_SEC = nothing
const SCIP_VERBLEVEL = 4

const MAX_MODEL_TRIES = 100


function usage()
    println("Usage: julia --project=. scripts/solve_random_qp_scip.jl [options]")
    println()
    println("Options:")
    println("  --parity             Run parity presolve before SCIP")
    println("  --no-parity          Solve the generated model directly with SCIP")
    println("  --nvars N            Number of variables (default: $N_VARS)")
    println("  --nconstraints M     Number of constraints (default: $N_CONSTRAINTS)")
    println("  --seed SEED          Random seed (default: $SEED)")
    println("  --density D          Quadratic term density (default: $DENSITY)")
    println("  --var-lb LB          Variable lower bound (default: $VAR_LB)")
    println("  --var-ub UB          Variable upper bound (default: $VAR_UB)")
    println("  --time-limit S       SCIP time limit in seconds, or none (default: unlimited)")
    println("  -h, --help           Show this help")
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
        nvars = N_VARS,
        nconstraints = N_CONSTRAINTS,
        var_lb = VAR_LB,
        var_ub = VAR_UB,
        seed = SEED,
        density = DENSITY,
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
        elseif arg == "--nvars"
            config = merge(config, (; nvars = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--nconstraints"
            config = merge(config, (; nconstraints = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--seed"
            config = merge(config, (; seed = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--density"
            config = merge(config, (; density = parse_value(Float64, args, idx, arg)))
            idx += 1
        elseif arg == "--var-lb"
            config = merge(config, (; var_lb = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--var-ub"
            config = merge(config, (; var_ub = parse_value(Int, args, idx, arg)))
            idx += 1
        elseif arg == "--time-limit" || arg == "--time-limit-sec"
            config = merge(config, (; time_limit_sec = parse_optional_float(args, idx, arg)))
            idx += 1
        else
            error("Unknown argument: $arg")
        end

        idx += 1
    end

    validate_config(config)
    return config
end


function validate_config(config)
    config.nvars > 0 || error("--nvars must be positive")
    config.nconstraints > 0 || error("--nconstraints must be positive")
    config.var_lb <= config.var_ub || error("--var-lb must be <= --var-ub")
    0.0 <= config.density <= 1.0 || error("--density must be in [0, 1]")
    return nothing
end


function random_quadratic_terms(
    rng::AbstractRNG,
    nvars::Int,
    density::Float64;
    force_odd_diagonal::Bool = false,
)::Vector{QuadTerm}
    terms = QuadTerm[]

    for i in 1:nvars
        rand(rng) < density || continue
        coeff = rand(rng, -500:500)
        coeff == 0 && continue
        push!(terms, (Float64(coeff), i, i))
    end

    for i in 1:(nvars - 1)
        for j in (i + 1):nvars
            rand(rng) < density || continue
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
        rng,
        nvars,
        density;
        force_odd_diagonal = force_odd_diagonal,
    )
    return terms
end


function build_random_constraint(
    rng::AbstractRNG,
    con_id::Int,
    x_star::Vector{Int},
    density::Float64;
    force_odd_diagonal::Bool = false,
)::PC.Constraint
    quad_terms = random_quadratic_terms(
        rng,
        length(x_star),
        density;
        force_odd_diagonal = force_odd_diagonal,
    )
    qe = PC.QuadExpr(quad_terms, LinTerm[]; constant = 0.0)
    rhs = Float64(PC.eval_full(qe, x_star))
    return PC.Constraint(con_id, qe, rhs, rhs)
end


function build_random_qp_model(config)
    rng = MersenneTwister(config.seed)
    x_star = rand(rng, config.var_lb:config.var_ub, config.nvars)

    vars = Dict(
        i => PC.IntVar(Float64(config.var_lb), Float64(config.var_ub))
        for i in 1:config.nvars
    )
    obj = PC.QuadExpr(QuadTerm[], LinTerm[]; constant = 0.0)

    for _ in 1:MAX_MODEL_TRIES
        cons = Vector{PC.Constraint}(undef, config.nconstraints)
        for con_id in 1:config.nconstraints
            cons[con_id] = build_random_constraint(
                rng,
                con_id,
                x_star,
                config.density;
                force_odd_diagonal = con_id == 1,
            )
        end

        model = PC.QPModel(copy(vars), cons, obj, :min)
        parity_model = PC.build_parity_model(model)
        if !isempty(parity_model.pos_to_var_id)
            return model
        end
    end

    error(
        "Failed to generate a model with nontrivial parity structure after $MAX_MODEL_TRIES attempts.",
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

    return (
        nvars = length(model.vars),
        total_constraints = total_constraints,
        quadratic_constraints = quadratic_constraints,
        linear_constraints = linear_constraints,
        infeasible = model.infeasible,
    )
end


function print_model_stats(label::AbstractString, stats)
    println(
        "$label: vars=$(stats.nvars) total_constraints=$(stats.total_constraints) " *
        "quadratic_constraints=$(stats.quadratic_constraints) " *
        "linear_constraints=$(stats.linear_constraints) infeasible=$(stats.infeasible)",
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
    qp_model = build_random_qp_model(config)
    original_stats = model_stats(qp_model)

    println(
        "nvars=$(config.nvars) nconstraints=$(config.nconstraints) " *
        "seed=$(config.seed) density=$(config.density) " *
        "var_lb=$(config.var_lb) var_ub=$(config.var_ub)",
    )
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
