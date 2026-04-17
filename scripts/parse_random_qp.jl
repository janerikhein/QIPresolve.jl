using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random

import QIPresolve as QIP

PC = QIP.PresolvingCore

const QuadTerm = Tuple{Float64, Int, Int}
const LinTerm = Tuple{Float64, Int}

const N_VARS = 30
const N_CONS = 30
const VAR_LB = 0
const VAR_UB = 100
const SEED = 17
const MAX_MODEL_TRIES = 100


function random_quadratic_terms(
        rng::AbstractRNG, nvars::Int; force_odd_diagonal::Bool = false
    )::Vector{QuadTerm}
    terms = QuadTerm[]

    for i in 1:nvars
        rand(rng) < 0.03 || continue
        coeff = rand(rng, -5:5)
        coeff == 0 && continue
        push!(terms, (Float64(coeff), i, i))
    end

    for i in 1:(nvars - 1)
        for j in (i + 1):nvars
            rand(rng) < 0.03 || continue
            coeff = rand(rng, -5:5)*2
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

    isempty(terms) && return random_quadratic_terms(rng, nvars; force_odd_diagonal = force_odd_diagonal)
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


function build_random_qp_model(rng::AbstractRNG)
    x_star = rand(rng, VAR_LB:VAR_UB, N_VARS)

    vars = Dict(i => PC.IntVar(Float64(VAR_LB), Float64(VAR_UB)) for i in 1:N_VARS)
    obj = PC.QuadExpr(QuadTerm[], LinTerm[]; constant = 0.0)

    for _ in 1:MAX_MODEL_TRIES
        cons = Vector{PC.Constraint}(undef, N_CONS)
        for con_id in 1:N_CONS
            cons[con_id] = build_random_constraint(
                rng, con_id, x_star; force_odd_diagonal = con_id == 1
            )
        end

        model = PC.QPModel(copy(vars), cons, obj, :min)
        xor_model = PC.build_parity_model(model)
        if !isempty(xor_model.pos_to_var_id)
            return model, x_star
        end
    end

    error("Failed to generate a model with nontrivial parity structure after $MAX_MODEL_TRIES attempts.")
end


function validate_model(qp_model::PC.QPModel, x_star::Vector{Int})
    @assert all(con -> con.lhs == con.rhs, qp_model.cons)
    @assert all(var -> var.lb == Float64(VAR_LB) && var.ub == Float64(VAR_UB), values(qp_model.vars))

    for con in qp_model.cons
        value = PC.eval_full(con.qe, x_star)
        @assert value == con.rhs
    end

    return nothing
end

solution_parity(x_star::Vector{Int}, vid::Int) = isodd(x_star[vid])

parity_label(parity::Bool) = parity ? "odd" : "even"

function collect_fixed_parity_matches(
        xor_model::PC.ParityModel,
        prop::PC.PropagationManager,
        x_star::Vector{Int},
    )
    matches = NamedTuple{(:vid, :propagated, :x_star, :match), Tuple{Int, Bool, Bool, Bool}}[]

    for vid in xor_model.pos_to_var_id
        propagated = PC.fixed_value(prop, vid)
        propagated === nothing && continue

        x_star_parity = solution_parity(x_star, vid)
        push!(matches, (
            vid = vid,
            propagated = propagated,
            x_star = x_star_parity,
            match = propagated == x_star_parity,
        ))
    end

    return matches
end

log_domain_sum(model::PC.QPModel) = sum(log(var.ub - var.lb) for var in values(model.vars); init = 0.0)

function build_reduced_solution(model::PC.QPModel)
    solution = Dict{Int, Float64}()

    for (var_id, var) in model.vars
        @assert var.lb == var.ub
        solution[var_id] = var.lb
    end

    return solution
end

function dense_solution(solution::Dict{Int, Float64}, nvars::Int)
    return [round(Int, solution[var_id]) for var_id in 1:nvars]
end

function satisfies_model(model::PC.QPModel, x::Vector{Int})
    return all(con -> PC.eval_full(con.qe, x) == con.rhs, model.cons)
end

function mismatches(x_ref::Vector{Int}, x_candidate::Vector{Int}; limit::Int = 10)
    diffs = NamedTuple{(:vid, :x_star, :reconstructed), Tuple{Int, Int, Int}}[]

    for vid in eachindex(x_ref)
        x_ref[vid] == x_candidate[vid] && continue
        push!(diffs, (vid = vid, x_star = x_ref[vid], reconstructed = x_candidate[vid]))
        length(diffs) >= limit && break
    end

    return diffs
end


rng = MersenneTwister(SEED)
qp_model, x_star = build_random_qp_model(rng)
validate_model(qp_model, x_star)
original_model = deepcopy(qp_model)
postsolver = PC.ParityPostsolver(keys(qp_model.vars))

println("x_star = $(collect(x_star))")
println("nvars=$(length(qp_model.vars)) ncons=$(length(qp_model.cons)) seed=$SEED")
println("log_domain_sum_before = $(log_domain_sum(qp_model))")
propagator = PC.PropagationManager(Int[])
phase_idx = 0

#println(qp_model)

while !qp_model.infeasible
    
    global phase_idx += 1
    stats = PC.parity_presolve_phase!(qp_model, propagator, postsolver)
    PC.scale_constraints_gcd!(qp_model)
    #println(qp_model)
    println(
        "phase $phase_idx: changed=$(stats.changed) " *
        "fixed_parities=$(stats.fixed_parities) " *
        "pattern_rewritten_vars=$(stats.pattern_rewritten_vars) " *
        "infeasible=$(qp_model.infeasible) " *
        "nvars=$(length(qp_model.vars)) ncons=$(length(qp_model.cons))"
    )

    stats.changed || break
end

#println(qp_model)

println("log_domain_sum_after = $(log_domain_sum(qp_model))")
println("final model infeasible = $(qp_model.infeasible)")

if qp_model.infeasible
    println("postsolve skipped: reduced model is infeasible")
else
    reduced_solution = build_reduced_solution(qp_model)
    reconstructed_solution = PC.postsolve(postsolver, reduced_solution)
    x_reconstructed = dense_solution(reconstructed_solution, length(x_star))

    matches_x_star = x_reconstructed == collect(x_star)
    reconstructed_satisfies_original = satisfies_model(original_model, x_reconstructed)
    mismatch_rows = mismatches(collect(x_star), x_reconstructed)

    println("reduced_solution = $reduced_solution")
    println("x_reconstructed = $x_reconstructed")
    println("matches_x_star = $matches_x_star")
    println("reconstructed_satisfies_original = $reconstructed_satisfies_original")
    println("nmismatched = $(count(!iszero, x_reconstructed .- collect(x_star)))")
    println("mismatch_examples = $mismatch_rows")
end

