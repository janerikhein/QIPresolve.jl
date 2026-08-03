using Random
using JuMP

const _RANDOM_QIP_MAX_TERM_TRIES = 1_000

"""
    generate_random_qip_model(nvars, ncons; kwargs...) -> (model, x_star)

Generate a randomized integer quadratic JuMP model and the sampled reference
point `x_star`. When `force_feasibility=true`, `x_star` is guaranteed to
satisfy all variable and constraint bounds.
"""
function generate_random_qip_model(
        nvars::Int,
        ncons::Int;
        p_con_eq::Real,
        var_threshold_lb::Integer,
        var_threshold_ub::Integer,
        p_var_is_candidate::Real,
        p_var_bilin::Real,
        p_var_diag::Real,
        p_var_lin::Real,
        coeff_lb::Integer,
        coeff_ub::Integer,
        force_diag_even::Bool,
        force_lin_even::Bool,
        force_feasibility::Bool,
        constraint_slack_range::AbstractVector{<:Integer},
        seed::Int = 0,
    )
    config = _validate_random_qip_config(
        nvars,
        ncons,
        p_con_eq,
        var_threshold_lb,
        var_threshold_ub,
        p_var_is_candidate,
        p_var_bilin,
        p_var_diag,
        p_var_lin,
        coeff_lb,
        coeff_ub,
        force_diag_even,
        force_lin_even,
        force_feasibility,
        constraint_slack_range,
    )

    rng = seed == 0 ? Random.default_rng() : MersenneTwister(seed)
    domain = config.var_threshold_lb:config.var_threshold_ub
    x_star = [rand(rng, domain) for _ in 1:nvars]
    var_lbs, var_ubs = _sample_random_qip_variable_bounds(
        rng,
        x_star,
        domain,
        config.force_feasibility,
    )

    model = Model()
    @variable(
        model,
        x[i in 1:nvars],
        Int,
        lower_bound = Float64(var_lbs[i]),
        upper_bound = Float64(var_ubs[i]),
    )

    objective_terms = _sample_random_qip_objective_terms(rng, nvars, config)
    if isempty(objective_terms)
        @objective(model, Min, 0)
    else
        objective_expr = @expression(
            model,
            sum(coeff * x[var_id] for (coeff, var_id) in objective_terms),
        )
        @objective(model, Min, objective_expr)
    end

    for _ in 1:ncons
        quad_terms, lin_terms = _sample_random_qip_constraint_terms(rng, nvars, config)
        reference_value = _eval_random_qip_terms(quad_terms, lin_terms, x_star)
        lhs, rhs = _sample_random_qip_constraint_bounds(rng, reference_value, config)
        constraint_expr = @expression(
            model,
            sum(coeff * x[var_id_1] * x[var_id_2] for (coeff, var_id_1, var_id_2) in quad_terms) +
            sum(coeff * x[var_id] for (coeff, var_id) in lin_terms),
        )

        if lhs == rhs
            @constraint(model, constraint_expr == Float64(lhs))
        else
            @constraint(model, Float64(lhs) <= constraint_expr <= Float64(rhs))
        end
    end

    return model, x_star
end

function _validate_random_qip_config(
        nvars::Int,
        ncons::Int,
        p_con_eq::Real,
        var_threshold_lb::Integer,
        var_threshold_ub::Integer,
        p_var_is_candidate::Real,
        p_var_bilin::Real,
        p_var_diag::Real,
        p_var_lin::Real,
        coeff_lb::Integer,
        coeff_ub::Integer,
        force_diag_even::Bool,
        force_lin_even::Bool,
        force_feasibility::Bool,
        constraint_slack_range::AbstractVector{<:Integer},
    )
    nvars > 0 || throw(ArgumentError("nvars must be positive, got $nvars"))
    ncons >= 0 || throw(ArgumentError("ncons must be nonnegative, got $ncons"))
    _validate_probability(:p_con_eq, p_con_eq)
    _validate_probability(:p_var_is_candidate, p_var_is_candidate)
    _validate_probability(:p_var_bilin, p_var_bilin)
    _validate_probability(:p_var_diag, p_var_diag)
    _validate_probability(:p_var_lin, p_var_lin)

    var_threshold_lb <= var_threshold_ub ||
        throw(ArgumentError("var_threshold_lb must be <= var_threshold_ub"))
    coeff_lb <= coeff_ub || throw(ArgumentError("coeff_lb must be <= coeff_ub"))
    _has_nonzero_coefficient(coeff_lb, coeff_ub) ||
        throw(ArgumentError("coefficient range must contain a nonzero value"))
    if force_diag_even && p_var_diag > 0
        _has_even_nonzero_coefficient(coeff_lb, coeff_ub) ||
            throw(ArgumentError("coefficient range must contain a nonzero even value when force_diag_even=true"))
    end
    if force_lin_even && p_var_lin > 0
        _has_even_nonzero_coefficient(coeff_lb, coeff_ub) ||
            throw(ArgumentError("coefficient range must contain a nonzero even value when force_lin_even=true"))
    end

    slack_values = Int.(collect(constraint_slack_range))
    isempty(slack_values) && throw(ArgumentError("constraint_slack_range must be nonempty"))
    if force_feasibility && p_con_eq < 1.0
        any(<=(0), slack_values) ||
            throw(ArgumentError("constraint_slack_range must contain a nonpositive value for feasible interval constraints"))
        any(>=(0), slack_values) ||
            throw(ArgumentError("constraint_slack_range must contain a nonnegative value for feasible interval constraints"))
    end

    return (
        p_con_eq = Float64(p_con_eq),
        var_threshold_lb = Int(var_threshold_lb),
        var_threshold_ub = Int(var_threshold_ub),
        p_var_is_candidate = Float64(p_var_is_candidate),
        p_var_bilin = Float64(p_var_bilin),
        p_var_diag = Float64(p_var_diag),
        p_var_lin = Float64(p_var_lin),
        coeff_lb = Int(coeff_lb),
        coeff_ub = Int(coeff_ub),
        force_diag_even = force_diag_even,
        force_lin_even = force_lin_even,
        force_feasibility = force_feasibility,
        slack_values = slack_values,
        feasible_lower_offsets = filter(<=(0), slack_values),
        feasible_upper_offsets = filter(>=(0), slack_values),
    )
end

function _validate_probability(name::Symbol, value::Real)
    value_float = Float64(value)
    isfinite(value_float) && 0.0 <= value_float <= 1.0 ||
        throw(ArgumentError("$name must be in [0, 1], got $value"))
    return nothing
end

function _has_nonzero_coefficient(coeff_lb::Integer, coeff_ub::Integer)
    return coeff_lb < 0 || coeff_ub > 0
end

function _has_even_nonzero_coefficient(coeff_lb::Integer, coeff_ub::Integer)
    for value in coeff_lb:coeff_ub
        value != 0 && iseven(value) && return true
    end
    return false
end

function _sample_random_qip_variable_bounds(
        rng::AbstractRNG,
        x_star::Vector{Int},
        domain::AbstractUnitRange{Int},
        force_feasibility::Bool,
    )
    lbs = Vector{Int}(undef, length(x_star))
    ubs = Vector{Int}(undef, length(x_star))

    for (idx, value) in enumerate(x_star)
        if force_feasibility
            lbs[idx] = rand(rng, first(domain):value)
            ubs[idx] = rand(rng, value:last(domain))
        else
            first_endpoint = rand(rng, domain)
            second_endpoint = rand(rng, domain)
            lbs[idx] = min(first_endpoint, second_endpoint)
            ubs[idx] = max(first_endpoint, second_endpoint)
        end
    end

    return lbs, ubs
end

function _sample_random_qip_objective_terms(rng::AbstractRNG, nvars::Int, config)
    terms = Tuple{Int, Int}[]
    sizehint!(terms, nvars)

    for var_id in 1:nvars
        rand(rng) < config.p_var_lin || continue
        push!(
            terms,
            (
                _sample_random_qip_coefficient(
                    rng,
                    config.coeff_lb,
                    config.coeff_ub;
                    force_even = config.force_lin_even,
                ),
                var_id,
            ),
        )
    end

    if isempty(terms) && config.p_var_lin > 0.0
        var_id = rand(rng, 1:nvars)
        push!(
            terms,
            (
                _sample_random_qip_coefficient(
                    rng,
                    config.coeff_lb,
                    config.coeff_ub;
                    force_even = config.force_lin_even,
                ),
                var_id,
            ),
        )
    end

    return terms
end

function _sample_random_qip_constraint_terms(rng::AbstractRNG, nvars::Int, config)
    for _ in 1:_RANDOM_QIP_MAX_TERM_TRIES
        candidates = _sample_random_qip_candidates(rng, nvars, config.p_var_is_candidate)
        quad_terms = Tuple{Int, Int, Int}[]
        lin_terms = Tuple{Int, Int}[]

        for var_id in candidates
            rand(rng) < config.p_var_diag || continue
            push!(
                quad_terms,
                (
                    _sample_random_qip_coefficient(
                        rng,
                        config.coeff_lb,
                        config.coeff_ub;
                        force_even = config.force_diag_even,
                    ),
                    var_id,
                    var_id,
                ),
            )
        end

        for (idx, first_id) in enumerate(candidates)
            for second_id in @view candidates[(idx + 1):end]
                rand(rng) < config.p_var_bilin || continue
                push!(
                    quad_terms,
                    (
                        _sample_random_qip_coefficient(
                            rng,
                            config.coeff_lb,
                            config.coeff_ub;
                            force_even = false,
                        ),
                        first_id,
                        second_id,
                    ),
                )
            end
        end

        for var_id in candidates
            rand(rng) < config.p_var_lin || continue
            push!(
                lin_terms,
                (
                    _sample_random_qip_coefficient(
                        rng,
                        config.coeff_lb,
                        config.coeff_ub;
                        force_even = config.force_lin_even,
                    ),
                    var_id,
                ),
            )
        end

        (!isempty(quad_terms) || !isempty(lin_terms)) && return quad_terms, lin_terms
    end

    throw(ArgumentError("failed to sample a nonempty constraint after $_RANDOM_QIP_MAX_TERM_TRIES attempts"))
end

function _sample_random_qip_candidates(
        rng::AbstractRNG,
        nvars::Int,
        p_var_is_candidate::Float64,
    )
    candidates = Int[]
    sizehint!(candidates, nvars)

    for var_id in 1:nvars
        rand(rng) < p_var_is_candidate || continue
        push!(candidates, var_id)
    end

    isempty(candidates) && push!(candidates, rand(rng, 1:nvars))
    return candidates
end

function _sample_random_qip_coefficient(
        rng::AbstractRNG,
        coeff_lb::Int,
        coeff_ub::Int;
        force_even::Bool,
    )
    for _ in 1:100
        coeff = rand(rng, coeff_lb:coeff_ub)
        coeff == 0 && continue
        force_even && isodd(coeff) && continue
        return coeff
    end

    valid_coefficients = [
        coeff for coeff in coeff_lb:coeff_ub
        if coeff != 0 && (!force_even || iseven(coeff))
    ]
    isempty(valid_coefficients) && throw(ArgumentError("coefficient range has no usable coefficient"))
    return rand(rng, valid_coefficients)
end

function _eval_random_qip_terms(
        quad_terms::Vector{Tuple{Int, Int, Int}},
        lin_terms::Vector{Tuple{Int, Int}},
        x_star::Vector{Int},
    )
    value = 0
    for (coeff, var_id_1, var_id_2) in quad_terms
        value += coeff * x_star[var_id_1] * x_star[var_id_2]
    end
    for (coeff, var_id) in lin_terms
        value += coeff * x_star[var_id]
    end
    return value
end

function _sample_random_qip_constraint_bounds(
        rng::AbstractRNG,
        reference_value::Int,
        config,
    )
    if rand(rng) < config.p_con_eq
        offset = config.force_feasibility ? 0 : rand(rng, config.slack_values)
        value = reference_value + offset
        return value, value
    elseif config.force_feasibility
        lower_offset = rand(rng, config.feasible_lower_offsets)
        upper_offset = rand(rng, config.feasible_upper_offsets)
        return reference_value + lower_offset, reference_value + upper_offset
    else
        first_offset = rand(rng, config.slack_values)
        second_offset = rand(rng, config.slack_values)
        lhs = reference_value + min(first_offset, second_offset)
        rhs = reference_value + max(first_offset, second_offset)
        return lhs, rhs
    end
end
