module ResidueConstraintExperimentCommon

using CSV
using Printf
using Random

import QIPresolve as QIP
import QIPresolve.PresolvingCore as PC

const DEFAULT_TREEWIDTH_THRESHOLD = QIP.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD
const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}

struct StrategySpec
    name::String
    moduli::Vector{Int}
end

const MODULI_FAMILY_STRATEGY_NAMES = [
    "primes_lt_64",
    "prime_powers_lt_64",
    "all_2_to_63",
    "primes_lt_128",
    "prime_powers_lt_128",
    "all_2_to_128",
]

struct ExactBoundTightening
    lhs::Float64
    rhs::Float64
    relative_bound_reduction::Float64
    assignment_count::Int
end

struct ExpressionTerms
    lin_terms::Vector{Tuple{PC.VarId, Int}}
    quad_terms::Vector{Tuple{PC.VarId, PC.VarId, Int}}
end

function parse_int(value::AbstractString, name::AbstractString)::Int
    parsed = tryparse(Int, strip(value))
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function parse_float(value::AbstractString, name::AbstractString)::Float64
    parsed = tryparse(Float64, strip(value))
    parsed === nothing && error("Invalid $name: $value")
    return parsed
end

function parse_bool(value::AbstractString, name::AbstractString)::Bool
    normalized = lowercase(strip(value))
    normalized in ("true", "1", "yes", "y") && return true
    normalized in ("false", "0", "no", "n") && return false
    error("Invalid $name: $value. Expected true/false, 1/0, or yes/no.")
end

function parse_int_list(value::AbstractString, name::AbstractString)::Vector{Int}
    normalized = replace(strip(value), " " => "")
    values = Int[]

    for part in split(normalized, ","; keepempty = false)
        if occursin(":", part)
            endpoints = split(part, ":"; limit = 2)
            length(endpoints) == 2 || error("Invalid $name range: $part")
            lower = parse_int(endpoints[1], name)
            upper = parse_int(endpoints[2], name)
            lower <= upper || error("Invalid $name range: lower bound must be <= upper bound")
            append!(values, lower:upper)
        else
            push!(values, parse_int(part, name))
        end
    end

    isempty(values) && error("$name must contain at least one integer")
    return values
end

function parse_string_list(value::AbstractString, name::AbstractString)::Vector{String}
    normalized = replace(strip(value), " " => "")
    values = String[]
    for part in split(normalized, ","; keepempty = false)
        push!(values, String(part))
    end

    isempty(values) && error("$name must contain at least one value")
    return values
end

function validate_probability(name::AbstractString, value::Float64)
    isfinite(value) && 0.0 <= value <= 1.0 ||
        error("$name must be in [0, 1], got $value")
    return value
end

function validate_nonnegative_alpha(alpha::Real)::Float64
    alpha_float = Float64(alpha)
    isfinite(alpha_float) || error("alpha must be finite, got $alpha")
    alpha_float >= 0.0 || error("alpha must be >= 0, got $alpha")
    return alpha_float
end

function coefficient_values(coeff_lb::Int, coeff_ub::Int)
    values = [value for value in coeff_lb:coeff_ub if value != 0]
    isempty(values) && error("coefficient range must contain a nonzero value")
    return values
end

function coefficient_palette(
        rng::AbstractRNG,
        coeff_lb::Int,
        coeff_ub::Int,
        max_distinct_coeffs::Int,
    )
    coefficients = coefficient_values(coeff_lb, coeff_ub)
    palette_size = min(max_distinct_coeffs, length(coefficients))
    palette_size == length(coefficients) && return coefficients

    order = randperm(rng, length(coefficients))
    palette = coefficients[order[1:palette_size]]
    return sort!(palette)
end

function is_prime(candidate::Int)
    candidate < 2 && return false
    candidate == 2 && return true
    iseven(candidate) && return false

    divisor = 3
    while divisor <= div(candidate, divisor)
        candidate % divisor == 0 && return false
        divisor += 2
    end
    return true
end

function primes_less_than(limit::Int)
    return [candidate for candidate in 2:(limit - 1) if is_prime(candidate)]
end

function prime_power_moduli_less_than(limit::Int)
    values = Set{Int}()
    for prime in primes_less_than(limit)
        value = prime
        while value < limit
            push!(values, value)
            value > div(limit - 1, prime) && break
            value *= prime
        end
    end
    return sort!(collect(values))
end

function full_moduli_less_than(limit::Int)
    upper = limit - 1
    return 2 <= upper ? collect(2:upper) : Int[]
end

function full_moduli_through(upper::Int)
    return 2 <= upper ? collect(2:upper) : Int[]
end

moduli_family_strategy_names() = copy(MODULI_FAMILY_STRATEGY_NAMES)

function normalize_moduli_strategy(value::AbstractString)::String
    normalized = lowercase(replace(strip(value), "-" => "_"))
    normalized == "primes_lt_64" && return normalized
    normalized == "primes_lt_128" && return normalized
    normalized == "prime_powers_lt_64" && return normalized
    normalized == "prime_powers_lt_128" && return normalized
    normalized == "all_2_to_63" && return normalized
    normalized == "all_2_to_128" && return normalized
    error(
        "Invalid moduli_strategy: $value. Expected one of " *
        join(MODULI_FAMILY_STRATEGY_NAMES, ", ") * ".",
    )
end

function normalize_moduli_strategies(values::AbstractVector{<:AbstractString})
    strategies = [normalize_moduli_strategy(value) for value in values]
    isempty(strategies) && error("moduli_strategies must contain at least one value")
    return strategies
end

function strategy_spec(name::AbstractString)
    normalized = normalize_moduli_strategy(name)
    normalized == "primes_lt_64" &&
        return StrategySpec(normalized, primes_less_than(64))
    normalized == "prime_powers_lt_64" &&
        return StrategySpec(normalized, prime_power_moduli_less_than(64))
    normalized == "all_2_to_63" &&
        return StrategySpec(normalized, full_moduli_less_than(64))
    normalized == "primes_lt_128" &&
        return StrategySpec(normalized, primes_less_than(128))
    normalized == "prime_powers_lt_128" &&
        return StrategySpec(normalized, prime_power_moduli_less_than(128))
    normalized == "all_2_to_128" &&
        return StrategySpec(normalized, full_moduli_through(128))
    error("unsupported moduli strategy: $name")
end

function strategy_specs(names::AbstractVector{<:AbstractString} = ["primes_lt_128"])
    return [strategy_spec(name) for name in normalize_moduli_strategies(names)]
end

function build_one_constraint_model(
        nvars::Int,
        con::PC.Constraint,
        domain_lb::Int,
        domain_ub::Int,
    )
    vars = Dict{PC.VarId, PC.IntVar}()
    for var_id in 1:nvars
        vars[var_id] = PC.IntVar(Float64(domain_lb), Float64(domain_ub))
    end
    return PC.QPModel(vars, [con], PC.QuadExpr(QuadTerm[], LinTerm[]), :min)
end

function one_constraint_model(model::PC.QPModel, con::PC.Constraint)
    return PC.QPModel(model.vars, [con], model.obj_expr, model.obj_sense)
end

bound_snapshot(con::PC.Constraint) = (lhs = con.lhs, rhs = con.rhs)

function relative_bound_reduction(before, lhs::Real, rhs::Real)
    isfinite(before.lhs) && isfinite(before.rhs) || return 0.0

    baseline_range = before.rhs - before.lhs
    baseline_range > 0.0 || return 0.0

    improvement = 0.0
    isfinite(lhs) && (improvement += max(0.0, Float64(lhs) - before.lhs))
    isfinite(rhs) && (improvement += max(0.0, before.rhs - Float64(rhs)))
    return improvement / baseline_range
end

relative_bound_reduction(before, con::PC.Constraint) =
    relative_bound_reduction(before, con.lhs, con.rhs)

function relative_bound_range_reduction(before, con::PC.Constraint)
    isfinite(before.lhs) && isfinite(before.rhs) || return 0.0

    before_range = before.rhs - before.lhs
    before_range > 0.0 || return 0.0

    after_range = con.rhs - con.lhs
    return (before_range - after_range) / before_range
end

function tightened_bound_count(before, con::PC.Constraint)
    count = 0
    isfinite(before.lhs) && con.lhs > before.lhs && (count += 1)
    isfinite(before.rhs) && con.rhs < before.rhs && (count += 1)
    return count
end

constraint_tightened_to_equality(before, con::PC.Constraint) =
    before.lhs != before.rhs && con.lhs == con.rhs

function run_residue_strategy!(
        model::PC.QPModel,
        con::PC.Constraint,
        moduli::AbstractVector{<:Integer},
        treewidth_threshold::Int,
    )
    model.infeasible && return true

    status = PC._residue_constraint_status(con, model.vars)
    if status == :infeasible
        model.infeasible = true
        return true
    elseif status == :skip || isempty(moduli)
        return false
    end

    standardized = PC._standardize_residue_constraint(con, model.vars, model.obj_expr)
    cache = PC._ResidueCacheEntry[]

    for modulus_value in moduli
        modulus = Int(modulus_value)
        residue_result = PC._compute_achievable_residues(
            modulus,
            standardized.con,
            standardized.var_bounds;
            treewidth_threshold = treewidth_threshold,
        )
        residue_result.saturated && continue

        push!(cache, PC._ResidueCacheEntry(modulus, residue_result.residues))
        PC._reapply_residue_cache!(con, cache, standardized.constraint_shift)

        if con.lhs > con.rhs
            model.infeasible = true
            return true
        end
    end

    return false
end

function _integer_coefficient(value::Float64, description::String)
    isinteger(value) || throw(ArgumentError("$description must be integer-valued, got $value"))
    return trunc(Int, value)
end

function expression_terms(qe::PC.QuadExpr)
    var_ids = sort!(collect(PC.vars(qe)))
    lin_terms = Tuple{PC.VarId, Int}[]
    quad_terms = Tuple{PC.VarId, PC.VarId, Int}[]

    for (index, var_id) in enumerate(var_ids)
        lin_coeff = _integer_coefficient(
            PC.get_lin_coeff(qe, var_id),
            "linear coefficient for $var_id",
        )
        lin_coeff == 0 || push!(lin_terms, (var_id, lin_coeff))

        for other_id in @view var_ids[index:end]
            quad_coeff = _integer_coefficient(
                PC.get_quad_coeff(qe, var_id, other_id),
                "quadratic coefficient for ($var_id, $other_id)",
            )
            quad_coeff == 0 || push!(quad_terms, (var_id, other_id, quad_coeff))
        end
    end

    return ExpressionTerms(lin_terms, quad_terms)
end

function eval_terms(terms::ExpressionTerms, values::AbstractVector{Int})
    total = 0
    @inbounds for (var_id, coefficient) in terms.lin_terms
        total += coefficient * values[var_id]
    end
    @inbounds for (first_id, second_id, coefficient) in terms.quad_terms
        total += coefficient * values[first_id] * values[second_id]
    end
    return total
end

function exact_domains(var_bounds::Dict{PC.VarId, PC.IntVar})
    var_ids = sort!(collect(keys(var_bounds)))
    lbs = Int[]
    ubs = Int[]
    sizehint!(lbs, length(var_ids))
    sizehint!(ubs, length(var_ids))

    for var_id in var_ids
        var = var_bounds[var_id]
        isfinite(var.lb) || throw(ArgumentError("variable $var_id must have finite lower bound"))
        isfinite(var.ub) || throw(ArgumentError("variable $var_id must have finite upper bound"))
        isinteger(var.lb) || throw(ArgumentError("variable $var_id lower bound must be integer-valued"))
        isinteger(var.ub) || throw(ArgumentError("variable $var_id upper bound must be integer-valued"))
        var.lb <= var.ub || throw(ArgumentError("variable $var_id has inconsistent bounds"))
        push!(lbs, trunc(Int, var.lb))
        push!(ubs, trunc(Int, var.ub))
    end

    return var_ids, lbs, ubs
end

function assignment_count(lbs::AbstractVector{Int}, ubs::AbstractVector{Int})
    count = 1
    for (lb, ub) in zip(lbs, ubs)
        count *= ub - lb + 1
    end
    return count
end

function advance_assignment!(
        values::Vector{Int},
        var_ids::AbstractVector{PC.VarId},
        lbs::AbstractVector{Int},
        ubs::AbstractVector{Int},
    )
    for position in length(var_ids):-1:1
        var_id = var_ids[position]
        if values[var_id] < ubs[position]
            values[var_id] += 1
            return true
        end
        values[var_id] = lbs[position]
    end
    return false
end

function exact_bound_tightening(con::PC.Constraint, var_bounds::Dict{PC.VarId, PC.IntVar})
    isfinite(con.lhs) || throw(ArgumentError("constraint lower bound must be finite"))
    isfinite(con.rhs) || throw(ArgumentError("constraint upper bound must be finite"))
    con.lhs <= con.rhs || throw(ArgumentError("constraint has inconsistent bounds"))

    var_ids, lbs, ubs = exact_domains(var_bounds)
    total_assignments = assignment_count(lbs, ubs)
    terms = expression_terms(con.qe)
    values = zeros(Int, maximum(var_ids; init = 0))

    for (position, var_id) in enumerate(var_ids)
        values[var_id] = lbs[position]
    end

    lower_threshold = ceil(Int, con.lhs)
    upper_threshold = floor(Int, con.rhs)
    best_lhs = typemax(Int)
    best_rhs = typemin(Int)

    for _ in 1:total_assignments
        value = eval_terms(terms, values)
        if value >= lower_threshold && value < best_lhs
            best_lhs = value
        end
        if value <= upper_threshold && value > best_rhs
            best_rhs = value
        end
        advance_assignment!(values, var_ids, lbs, ubs)
    end

    best_lhs == typemax(Int) && error("no assignment attains a value >= lower bound")
    best_rhs == typemin(Int) && error("no assignment attains a value <= upper bound")

    before = bound_snapshot(con)
    lhs = Float64(best_lhs)
    rhs = Float64(best_rhs)
    return ExactBoundTightening(
        lhs,
        rhs,
        relative_bound_reduction(before, lhs, rhs),
        total_assignments,
    )
end

function bounds_fully_tightened_to_optimal_count(
        before,
        con::PC.Constraint,
        exact::ExactBoundTightening,
    )
    count = 0
    isfinite(before.lhs) && con.lhs > before.lhs && con.lhs == exact.lhs && (count += 1)
    isfinite(before.rhs) && con.rhs < before.rhs && con.rhs == exact.rhs && (count += 1)
    return count
end

function bound_gap_to_optimal(con::PC.Constraint, exact::ExactBoundTightening)
    return abs(con.lhs - exact.lhs) + abs(con.rhs - exact.rhs)
end

function _run_residue_trial!(
        result,
        model::PC.QPModel,
        con::PC.Constraint,
        treewidth_threshold::Int,
    )
    trial_con = deepcopy(con)
    trial_model = one_constraint_model(model, trial_con)
    before = bound_snapshot(trial_con)

    start_time = time()
    run_residue_strategy!(
        trial_model,
        trial_con,
        result.moduli,
        treewidth_threshold,
    )
    result.total_residue_time_sec += time() - start_time

    result.constraints += 1
    result.bounds_tightened += tightened_bound_count(before, trial_con)
    result.constraints_tightened_to_equality +=
        constraint_tightened_to_equality(before, trial_con) ? 1 : 0
    result.total_relative_bound_range_reduction +=
        relative_bound_range_reduction(before, trial_con)
    return trial_con, before
end

function record_strategy_trial!(
        result,
        model::PC.QPModel,
        con::PC.Constraint,
        exact::ExactBoundTightening,
        treewidth_threshold::Int,
    )
    trial_con, before = _run_residue_trial!(result, model, con, treewidth_threshold)
    result.bounds_fully_tightened_to_optimal +=
        bounds_fully_tightened_to_optimal_count(before, trial_con, exact)
    result.total_bound_gap_to_optimal += bound_gap_to_optimal(trial_con, exact)
    result.exact_assignments_per_constraint = exact.assignment_count

    return result
end

function record_strategy_trial!(
        result,
        model::PC.QPModel,
        con::PC.Constraint,
        ::Nothing,
        treewidth_threshold::Int,
    )
    _run_residue_trial!(result, model, con, treewidth_threshold)
    return result
end

rate(numerator::Real, denominator::Int) = denominator == 0 ? 0.0 : numerator / denominator

metric_float(value) = value === missing ? "missing" : @sprintf("%.6f", value)

function result_metric_tail(result)
    bounds_considered = 2 * result.constraints
    exact_enabled = getproperty(result, :exact_enumeration)
    return (
        strategy = result.name,
        moduli = join(result.moduli, " "),
        num_moduli = length(result.moduli),
        constraints = result.constraints,
        bounds_considered = bounds_considered,
        pct_constraints_tightened_to_equality = rate(
            100.0 * result.constraints_tightened_to_equality,
            result.constraints,
        ),
        avg_relative_bound_range_reduction = rate(
            result.total_relative_bound_range_reduction,
            result.constraints,
        ),
        exact_assignments_per_constraint = exact_enabled ?
            result.exact_assignments_per_constraint :
            missing,
        bounds_fully_tightened_to_optimal = exact_enabled ?
            result.bounds_fully_tightened_to_optimal :
            missing,
        pct_bounds_fully_tightened_to_optimal = exact_enabled ? rate(
            100.0 * result.bounds_fully_tightened_to_optimal,
            bounds_considered,
        ) : missing,
        bounds_tightened = result.bounds_tightened,
        pct_bounds_tightened = rate(100.0 * result.bounds_tightened, bounds_considered),
        avg_bound_gap_to_optimal = exact_enabled ? rate(
            result.total_bound_gap_to_optimal,
            bounds_considered,
        ) : missing,
        total_residue_time_sec = result.total_residue_time_sec,
        avg_wall_time_sec_per_constraint = rate(
            result.total_residue_time_sec,
            result.constraints,
        ),
    )
end

function write_csv(path::AbstractString, rows)
    mkpath(dirname(path))
    CSV.write(path, rows)
    return path
end

end # module
