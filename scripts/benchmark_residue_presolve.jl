#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Printf
using Random

import QIPresolve as QIP

PC = QIP.PresolvingCore

const QuadTerm = Tuple{Float64, Int, Int}
const LinTerm = Tuple{Float64, Int}

const N_VARS = 200
const N_CONSTRAINTS = 10
const N_INSTANCES = 10
const VAR_LB = -100
const VAR_UB = 100
const SEED = 31

const DIAG_DENSITY = 0.0
const BILINEAR_DENSITY = 0.00
const LINEAR_DENSITY = 0.01

const DIAG_COEFF_RANGE = -50:50
const BILINEAR_COEFF_RANGE = -50:50
const LINEAR_COEFF_RANGE = -50:50
const BOUND_SLACK_RANGE = 0:30

const MODULUS_STRATEGY = :divisor_free
const MODULUS_THRESHOLD = 256
const TREEWIDTH_THRESHOLD = 2
const RUN_WARMUP = true


function sample_nonzero(rng::AbstractRNG, values)
    coeff = rand(rng, values)
    while coeff == 0
        coeff = rand(rng, values)
    end
    return coeff
end


function random_terms(rng::AbstractRNG, nvars::Int)
    quad_terms = QuadTerm[]
    lin_terms = LinTerm[]

    for var_id in 1:nvars
        rand(rng) < DIAG_DENSITY || continue
        coeff = sample_nonzero(rng, DIAG_COEFF_RANGE)
        push!(quad_terms, (Float64(coeff), var_id, var_id))
    end

    for first_id in 1:(nvars - 1)
        for second_id in (first_id + 1):nvars
            rand(rng) < BILINEAR_DENSITY || continue
            coeff = 2 * sample_nonzero(rng, BILINEAR_COEFF_RANGE)
            push!(quad_terms, (Float64(coeff), first_id, second_id))
        end
    end

    for var_id in 1:nvars
        rand(rng) < LINEAR_DENSITY || continue
        coeff = sample_nonzero(rng, LINEAR_COEFF_RANGE)
        push!(lin_terms, (Float64(coeff), var_id))
    end

    if isempty(quad_terms) && isempty(lin_terms)
        var_id = rand(rng, 1:nvars)
        coeff = sample_nonzero(rng, LINEAR_COEFF_RANGE)
        push!(lin_terms, (Float64(coeff), var_id))
    end

    return quad_terms, lin_terms
end


function build_random_constraint(rng::AbstractRNG, con_id::Int, x_star::Vector{Int})
    quad_terms, lin_terms = random_terms(rng, length(x_star))
    qe = PC.QuadExpr(quad_terms, lin_terms; constant = 0.0)
    value = PC.eval_full(qe, x_star)
    isinteger(value) || error("generated non-integer feasible value $value")

    center = round(Int, value)
    lhs = center - rand(rng, BOUND_SLACK_RANGE)
    rhs = center + rand(rng, BOUND_SLACK_RANGE)

    return PC.Constraint(con_id, qe, Float64(lhs), Float64(rhs))
end


function build_random_model(rng::AbstractRNG)
    x_star = rand(rng, VAR_LB:VAR_UB, N_VARS)
    vars = Dict(
        var_id => PC.IntVar(Float64(VAR_LB), Float64(VAR_UB))
        for var_id in 1:N_VARS
    )
    cons = [
        build_random_constraint(rng, con_id, x_star)
        for con_id in 1:N_CONSTRAINTS
    ]
    obj = PC.QuadExpr(QuadTerm[], LinTerm[]; constant = 0.0)
    model = PC.QPModel(vars, cons, obj, :min)
    validate_feasible_point(model, x_star)
    #println(model)
    return model, x_star
end


function validate_feasible_point(model::PC.QPModel, x_star::Vector{Int})
    for con in model.cons
        value = PC.eval_full(con.qe, x_star)
        con.lhs <= value <= con.rhs ||
            error("sampled point violates constraint $(con.id): $(con.lhs) <= $value <= $(con.rhs)")
    end
    return nothing
end


function bound_improvement_stats(
        original::PC.QPModel,
        reduced::PC.QPModel,
    )
    improved_constraints = 0
    improved_bounds = 0
    total_improvement = 0.0

    for (old_con, new_con) in zip(original.cons, reduced.cons)
        constraint_improved = false

        lhs_improvement = new_con.lhs - old_con.lhs
        if lhs_improvement > 0.0
            improved_bounds += 1
            total_improvement += lhs_improvement
            constraint_improved = true
        end

        rhs_improvement = old_con.rhs - new_con.rhs
        if rhs_improvement > 0.0
            improved_bounds += 1
            total_improvement += rhs_improvement
            constraint_improved = true
        end

        constraint_improved && (improved_constraints += 1)
    end

    average_improvement = improved_bounds == 0 ? 0.0 : total_improvement / improved_bounds
    return (
        improved_constraints = improved_constraints,
        improved_bounds = improved_bounds,
        total_improvement = total_improvement,
        average_improvement = average_improvement,
    )
end


function print_config()
    println("Residue presolve benchmark")
    println("nvars = $N_VARS")
    println("nconstraints = $N_CONSTRAINTS")
    println("ninstances = $N_INSTANCES")
    println("var_domain = $VAR_LB:$VAR_UB")
    println("seed = $SEED")
    println("diag_density = $DIAG_DENSITY")
    println("bilinear_density = $BILINEAR_DENSITY")
    println("linear_density = $LINEAR_DENSITY")
    println("diag_coeff_range = $DIAG_COEFF_RANGE")
    println("bilinear_coeff_range = 2 * $BILINEAR_COEFF_RANGE")
    println("linear_coeff_range = $LINEAR_COEFF_RANGE")
    println("bound_slack_range = $BOUND_SLACK_RANGE")
    println("modulus_strategy = $MODULUS_STRATEGY")
    println("modulus_threshold = $MODULUS_THRESHOLD")
    println("treewidth_threshold = $TREEWIDTH_THRESHOLD")
    println("run_warmup = $RUN_WARMUP")
    println()
    return nothing
end


function run_warmup()
    rng = MersenneTwister(SEED - 1)
    model, _ = build_random_model(rng)
    PC.residue_presolve!(
        model,
        MODULUS_STRATEGY;
        threshold = MODULUS_THRESHOLD,
        treewidth_threshold = TREEWIDTH_THRESHOLD,
    )
    return nothing
end


function run_benchmark()
    rng = MersenneTwister(SEED)
    total_runtime = 0.0
    total_improved_constraints = 0
    total_improved_bounds = 0
    total_improvement = 0.0
    infeasible_count = 0

    print_config()
    if RUN_WARMUP
        println("running warmup instance")
        run_warmup()
        println()
    end
    println("instance\truntime_sec\timproved_constraints\timproved_bounds\tavg_bound_improvement\tinfeasible")

    for instance_idx in 1:N_INSTANCES
        original_model, x_star = build_random_model(rng)
        reduced_model = deepcopy(original_model)

        runtime = @elapsed PC.residue_presolve!(
            reduced_model,
            MODULUS_STRATEGY;
            threshold = MODULUS_THRESHOLD,
            treewidth_threshold = TREEWIDTH_THRESHOLD,
        )
        total_runtime += runtime

        reduced_model.infeasible && (infeasible_count += 1)
        reduced_model.infeasible || validate_feasible_point(reduced_model, x_star)

        stats = bound_improvement_stats(original_model, reduced_model)
        total_improved_constraints += stats.improved_constraints
        total_improved_bounds += stats.improved_bounds
        total_improvement += stats.total_improvement

        @printf(
            "%d\t%.6f\t%d\t%d\t%.6f\t%s\n",
            instance_idx,
            runtime,
            stats.improved_constraints,
            stats.improved_bounds,
            stats.average_improvement,
            string(reduced_model.infeasible),
        )
    end

    total_constraints = N_INSTANCES * N_CONSTRAINTS
    improved_fraction = total_constraints == 0 ? 0.0 : total_improved_constraints / total_constraints
    average_runtime = N_INSTANCES == 0 ? 0.0 : total_runtime / N_INSTANCES
    average_improvement = total_improved_bounds == 0 ? 0.0 : total_improvement / total_improved_bounds

    println()
    println("Aggregate")
    @printf("total_runtime_sec = %.6f\n", total_runtime)
    @printf("avg_runtime_sec_per_instance = %.6f\n", average_runtime)
    println("constraints_improved = $total_improved_constraints / $total_constraints")
    @printf("constraint_improvement_fraction = %.6f\n", improved_fraction)
    println("improved_individual_bounds = $total_improved_bounds")
    @printf("avg_bound_improvement_over_improved_bounds = %.6f\n", average_improvement)
    println("infeasible_count = $infeasible_count")

    return nothing
end


run_benchmark()
