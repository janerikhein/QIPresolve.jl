using Statistics: mean, median

Base.@kwdef struct ModelStats
    num_variables::Int
    num_binary_variables::Int
    num_general_integer_variables::Int
    num_fixed_variables::Int
    num_constraints::Int
    num_linear_constraints::Int
    num_linear_inequalities::Int
    num_linear_equalities::Int
    num_quadratic_constraints::Int
    num_quadratic_equalities::Int
    num_quadratic_inequalities::Int
    num_equality_constraints::Int
    num_inequality_constraints::Int
    num_linear_terms::Int
    num_quadratic_terms::Int
    num_diagonal_quadratic_terms::Int
    num_bilinear_terms::Int
    num_nonzeros::Int
    log_domain_sum::Float64
    max_domain_size::Int
    avg_domain_size::Float64
    median_domain_size::Float64
    num_distinct_coefficients::Int
    min_abs_coefficient::Float64
    max_abs_coefficient::Float64
    avg_abs_coefficient::Float64
    median_abs_coefficient::Float64
    coefficient_abs_range::Float64
    num_unit_coefficients::Int
    num_even_coefficients::Int
    num_odd_coefficients::Int
    num_coefficients_divisible_by_4::Int
    num_coefficients_2mod4::Int
    num_even_linear_coefficients::Float64
    num_even_diagonal_coefficients::Float64
    num_even_bilinear_coefficients::Float64
    num_divisible_by_4_linear_coefficients::Float64
    num_divisible_by_4_diagonal_coefficients::Float64
    num_constraints_all_linear_coefficients_even::Int
    num_constraints_all_diagonal_linear_coefficients_even::Int
    num_constraints_with_odd_diagonal_or_linear_coefficient::Int
    max_abs_linear_coefficient::Float64
    avg_abs_linear_coefficient::Float64
    max_abs_diagonal_coefficient::Float64
    avg_abs_diagonal_coefficient::Float64
    max_abs_bilinear_coefficient::Float64
    avg_abs_bilinear_coefficient::Float64
    min_constraint_bound_width::Float64
    max_constraint_bound_width::Float64
    avg_constraint_bound_width::Float64
    median_constraint_bound_width::Float64
    min_constraint_lower_bound::Float64
    max_constraint_lower_bound::Float64
    avg_constraint_lower_bound::Float64
    min_constraint_upper_bound::Float64
    max_constraint_upper_bound::Float64
    avg_constraint_upper_bound::Float64
    max_abs_constraint_bound::Float64
    avg_abs_constraint_bound::Float64
    median_abs_constraint_bound::Float64
    num_variables_with_nonfinite_domain::Int
    num_constraints_with_fractional_coefficients::Int
end

@inline _zero_if_empty_min(values::AbstractVector{<:Real}) = isempty(values) ? 0.0 : Float64(minimum(values))
@inline _zero_if_empty_max(values::AbstractVector{<:Real}) = isempty(values) ? 0.0 : Float64(maximum(values))
@inline _zero_if_empty_avg(values::AbstractVector{<:Real}) = isempty(values) ? 0.0 : Float64(mean(values))
@inline _zero_if_empty_median(values::AbstractVector{<:Real}) = isempty(values) ? 0.0 : Float64(median(values))

function _abs_values(values::AbstractVector{Float64})
    return Float64[abs(value) for value in values]
end

@inline _is_integer_coefficient(value::Float64) = isfinite(value) && isinteger(value)
@inline _is_even_coefficient(value::Float64) = _is_integer_coefficient(value) && mod(value, 2.0) == 0.0
@inline _is_odd_coefficient(value::Float64) = _is_integer_coefficient(value) && mod(value, 2.0) == 1.0
@inline _is_divisible_by_4_coefficient(value::Float64) = _is_integer_coefficient(value) && mod(value, 4.0) == 0.0
@inline _is_2mod4_coefficient(value::Float64) = _is_integer_coefficient(value) && mod(value, 4.0) == 2.0

function _nonzero_linear_coefficients(qe::QuadExpr)
    coeffs = Float64[]
    sizehint!(coeffs, qe.nvars)

    @inbounds for pos in 1:qe.nvars
        coeff = qe.lin_buf[qe.perm[pos]]
        coeff == 0.0 || push!(coeffs, coeff)
    end

    return coeffs
end

function _nonzero_quadratic_coefficients(qe::QuadExpr)
    diagonal_coeffs = Float64[]
    bilinear_coeffs = Float64[]
    sizehint!(diagonal_coeffs, qe.nvars)
    sizehint!(bilinear_coeffs, qe.nvars)

    @inbounds for pos_i in 1:qe.nvars
        phys_i = qe.perm[pos_i]
        diag_coeff = qe.quad_buf[phys_i, phys_i]
        diag_coeff == 0.0 || push!(diagonal_coeffs, diag_coeff)

        for pos_j in (pos_i + 1):qe.nvars
            phys_j = qe.perm[pos_j]
            coeff = if phys_i <= phys_j
                qe.quad_buf[phys_i, phys_j]
            else
                qe.quad_buf[phys_j, phys_i]
            end
            coeff == 0.0 || push!(bilinear_coeffs, coeff)
        end
    end

    return diagonal_coeffs, bilinear_coeffs
end

@inline function _has_fractional_coefficient(coeffs::AbstractVector{Float64})
    return any(coeff -> !_is_integer_coefficient(coeff), coeffs)
end

@inline function _all_even_coefficients(coeffs::AbstractVector{Float64})
    return all(_is_even_coefficient, coeffs)
end

function _finite_domain_sizes(model::QPModel)
    domain_sizes = Int[]
    sizehint!(domain_sizes, length(model.vars))
    nonfinite_domains = 0

    for var in values(model.vars)
        if !isfinite(var.lb) || !isfinite(var.ub)
            nonfinite_domains += 1
            continue
        end

        lower = ceil(Int, var.lb)
        upper = floor(Int, var.ub)
        push!(domain_sizes, max(0, upper - lower + 1))
    end

    return domain_sizes, nonfinite_domains
end

"""
    ModelStats(model)

Compute a non-mutating snapshot of structural statistics for `model`.

Only constraints in `model.cons` contribute to coefficient, term, and bound
statistics. Objective coefficients are intentionally excluded.
"""
function ModelStats(model::QPModel)
    domain_sizes, num_variables_with_nonfinite_domain = _finite_domain_sizes(model)
    log_domain_sum = sum((size > 0 ? log2(size) : 0.0 for size in domain_sizes); init = 0.0)

    linear_coefficients = Float64[]
    diagonal_coefficients = Float64[]
    bilinear_coefficients = Float64[]
    coefficients = Float64[]
    constraint_bound_widths = Float64[]
    constraint_lower_bounds = Float64[]
    constraint_upper_bounds = Float64[]
    abs_constraint_bounds = Float64[]

    num_linear_constraints = 0
    num_linear_equalities = 0
    num_linear_inequalities = 0
    num_quadratic_constraints = 0
    num_quadratic_equalities = 0
    num_quadratic_inequalities = 0
    num_equality_constraints = 0
    num_inequality_constraints = 0
    num_constraints_all_linear_coefficients_even = 0
    num_constraints_all_diagonal_linear_coefficients_even = 0
    num_constraints_with_odd_diagonal_or_linear_coefficient = 0
    num_constraints_with_fractional_coefficients = 0

    for con in model.cons
        con_linear_coeffs = _nonzero_linear_coefficients(con.qe)
        con_diagonal_coeffs, con_bilinear_coeffs = _nonzero_quadratic_coefficients(con.qe)
        con_quadratic_coeffs = [con_diagonal_coeffs; con_bilinear_coeffs]
        con_coeffs = [con_linear_coeffs; con_quadratic_coeffs]

        append!(linear_coefficients, con_linear_coeffs)
        append!(diagonal_coefficients, con_diagonal_coeffs)
        append!(bilinear_coefficients, con_bilinear_coeffs)
        append!(coefficients, con_coeffs)

        is_quadratic = !isempty(con_quadratic_coeffs)
        is_eq = is_equality(con)
        if is_quadratic
            num_quadratic_constraints += 1
            if is_eq
                num_quadratic_equalities += 1
            else
                num_quadratic_inequalities += 1
            end

            if _all_even_coefficients(con_linear_coeffs) && _all_even_coefficients(con_diagonal_coeffs)
                num_constraints_all_diagonal_linear_coefficients_even += 1
            end
            if any(_is_odd_coefficient, con_linear_coeffs) || any(_is_odd_coefficient, con_diagonal_coeffs)
                num_constraints_with_odd_diagonal_or_linear_coefficient += 1
            end
        else
            num_linear_constraints += 1
            if is_eq
                num_linear_equalities += 1
            else
                num_linear_inequalities += 1
            end
        end

        if is_eq
            num_equality_constraints += 1
        else
            num_inequality_constraints += 1
        end

        _all_even_coefficients(con_linear_coeffs) &&
            (num_constraints_all_linear_coefficients_even += 1)
        _has_fractional_coefficient(con_coeffs) &&
            (num_constraints_with_fractional_coefficients += 1)

        if isfinite(con.lhs)
            push!(constraint_lower_bounds, con.lhs)
            push!(abs_constraint_bounds, abs(con.lhs))
        end
        if isfinite(con.rhs)
            push!(constraint_upper_bounds, con.rhs)
            push!(abs_constraint_bounds, abs(con.rhs))
        end
        if isfinite(con.lhs) && isfinite(con.rhs)
            width = con.rhs - con.lhs
            width > 0.0 && push!(constraint_bound_widths, width)
        end
    end

    abs_coefficients = _abs_values(coefficients)
    abs_linear_coefficients = _abs_values(linear_coefficients)
    abs_diagonal_coefficients = _abs_values(diagonal_coefficients)
    abs_bilinear_coefficients = _abs_values(bilinear_coefficients)
    min_abs_coefficient = _zero_if_empty_min(abs_coefficients)
    max_abs_coefficient = _zero_if_empty_max(abs_coefficients)

    return ModelStats(
        num_variables = length(model.vars),
        num_binary_variables = count(is_binary, values(model.vars)),
        num_general_integer_variables = count(!is_binary(var) for var in values(model.vars)),
        num_fixed_variables = count(var -> var.lb == var.ub, values(model.vars)),
        num_constraints = length(model.cons),
        num_linear_constraints = num_linear_constraints,
        num_linear_inequalities = num_linear_inequalities,
        num_linear_equalities = num_linear_equalities,
        num_quadratic_constraints = num_quadratic_constraints,
        num_quadratic_equalities = num_quadratic_equalities,
        num_quadratic_inequalities = num_quadratic_inequalities,
        num_equality_constraints = num_equality_constraints,
        num_inequality_constraints = num_inequality_constraints,
        num_linear_terms = length(linear_coefficients),
        num_quadratic_terms = length(diagonal_coefficients) + length(bilinear_coefficients),
        num_diagonal_quadratic_terms = length(diagonal_coefficients),
        num_bilinear_terms = length(bilinear_coefficients),
        num_nonzeros = length(coefficients),
        log_domain_sum = Float64(log_domain_sum),
        max_domain_size = isempty(domain_sizes) ? 0 : maximum(domain_sizes),
        avg_domain_size = _zero_if_empty_avg(domain_sizes),
        median_domain_size = _zero_if_empty_median(domain_sizes),
        num_distinct_coefficients = length(Set(coefficients)),
        min_abs_coefficient = min_abs_coefficient,
        max_abs_coefficient = max_abs_coefficient,
        avg_abs_coefficient = _zero_if_empty_avg(abs_coefficients),
        median_abs_coefficient = _zero_if_empty_median(abs_coefficients),
        coefficient_abs_range = max_abs_coefficient - min_abs_coefficient,
        num_unit_coefficients = count(coeff -> abs(coeff) == 1.0, coefficients),
        num_even_coefficients = count(_is_even_coefficient, coefficients),
        num_odd_coefficients = count(_is_odd_coefficient, coefficients),
        num_coefficients_divisible_by_4 = count(_is_divisible_by_4_coefficient, coefficients),
        num_coefficients_2mod4 = count(_is_2mod4_coefficient, coefficients),
        num_even_linear_coefficients = Float64(count(_is_even_coefficient, linear_coefficients)),
        num_even_diagonal_coefficients = Float64(count(_is_even_coefficient, diagonal_coefficients)),
        num_even_bilinear_coefficients = Float64(count(_is_even_coefficient, bilinear_coefficients)),
        num_divisible_by_4_linear_coefficients = Float64(count(_is_divisible_by_4_coefficient, linear_coefficients)),
        num_divisible_by_4_diagonal_coefficients = Float64(count(_is_divisible_by_4_coefficient, diagonal_coefficients)),
        num_constraints_all_linear_coefficients_even = num_constraints_all_linear_coefficients_even,
        num_constraints_all_diagonal_linear_coefficients_even = num_constraints_all_diagonal_linear_coefficients_even,
        num_constraints_with_odd_diagonal_or_linear_coefficient = num_constraints_with_odd_diagonal_or_linear_coefficient,
        max_abs_linear_coefficient = _zero_if_empty_max(abs_linear_coefficients),
        avg_abs_linear_coefficient = _zero_if_empty_avg(abs_linear_coefficients),
        max_abs_diagonal_coefficient = _zero_if_empty_max(abs_diagonal_coefficients),
        avg_abs_diagonal_coefficient = _zero_if_empty_avg(abs_diagonal_coefficients),
        max_abs_bilinear_coefficient = _zero_if_empty_max(abs_bilinear_coefficients),
        avg_abs_bilinear_coefficient = _zero_if_empty_avg(abs_bilinear_coefficients),
        min_constraint_bound_width = _zero_if_empty_min(constraint_bound_widths),
        max_constraint_bound_width = _zero_if_empty_max(constraint_bound_widths),
        avg_constraint_bound_width = _zero_if_empty_avg(constraint_bound_widths),
        median_constraint_bound_width = _zero_if_empty_median(constraint_bound_widths),
        min_constraint_lower_bound = _zero_if_empty_min(constraint_lower_bounds),
        max_constraint_lower_bound = _zero_if_empty_max(constraint_lower_bounds),
        avg_constraint_lower_bound = _zero_if_empty_avg(constraint_lower_bounds),
        min_constraint_upper_bound = _zero_if_empty_min(constraint_upper_bounds),
        max_constraint_upper_bound = _zero_if_empty_max(constraint_upper_bounds),
        avg_constraint_upper_bound = _zero_if_empty_avg(constraint_upper_bounds),
        max_abs_constraint_bound = _zero_if_empty_max(abs_constraint_bounds),
        avg_abs_constraint_bound = _zero_if_empty_avg(abs_constraint_bounds),
        median_abs_constraint_bound = _zero_if_empty_median(abs_constraint_bounds),
        num_variables_with_nonfinite_domain = num_variables_with_nonfinite_domain,
        num_constraints_with_fractional_coefficients = num_constraints_with_fractional_coefficients,
    )
end
