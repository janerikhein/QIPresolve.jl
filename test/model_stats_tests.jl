using Test
import QIPresolve.PresolvingCore as PC

const StatsQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const StatsLinTerm = Tuple{Float64, PC.VarId}
const _MODEL_STATS_NEXT_CON_ID = Ref(0)

model_stats_next_con_id() = (_MODEL_STATS_NEXT_CON_ID[] += 1)
model_stats_empty_objective() = PC.QuadExpr(StatsQuadTerm[], StatsLinTerm[])

@testset "ModelStats computes structural model snapshot" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 3.0),
        3 => PC.IntVar(2.0, 2.0),
        4 => PC.IntVar(-Inf, 5.0),
        5 => PC.IntVar(1.2, 4.8),
    )

    linear_equality = PC.Constraint(
        model_stats_next_con_id(),
        PC.QuadExpr(StatsQuadTerm[], StatsLinTerm[(2.0, 1), (-4.0, 2)]),
        6.0,
        6.0,
    )
    linear_ranged = PC.Constraint(
        model_stats_next_con_id(),
        PC.QuadExpr(StatsQuadTerm[], StatsLinTerm[(3.0, 2), (5.0, 3)]),
        1.0,
        9.0,
    )
    quadratic_equality = PC.Constraint(
        model_stats_next_con_id(),
        PC.QuadExpr(StatsQuadTerm[(6.0, 1, 1), (8.0, 1, 2)], StatsLinTerm[(4.0, 2)]),
        10.0,
        10.0,
    )
    quadratic_ranged = PC.Constraint(
        model_stats_next_con_id(),
        PC.QuadExpr(StatsQuadTerm[(1.5, 2, 2), (2.0, 2, 3)], StatsLinTerm[(-3.0, 3)]),
        -2.0,
        8.0,
    )
    linear_one_sided = PC.Constraint(
        model_stats_next_con_id(),
        PC.QuadExpr(StatsQuadTerm[], StatsLinTerm[(7.0, 1)]),
        -Inf,
        10.0,
    )

    model = PC.QPModel(
        vars,
        [linear_equality, linear_ranged, quadratic_equality, quadratic_ranged, linear_one_sided],
        model_stats_empty_objective(),
        :min,
    )

    stats = PC.ModelStats(model)

    @test stats.num_variables == 5
    @test stats.num_binary_variables == 1
    @test stats.num_general_integer_variables == 4
    @test stats.num_fixed_variables == 1
    @test stats.num_variables_with_nonfinite_domain == 1
    @test stats.max_domain_size == 4
    @test stats.avg_domain_size == 2.5
    @test stats.median_domain_size == 2.5
    @test isapprox(stats.log_domain_sum, 3.0 + log2(3.0); atol = 1.0e-12)

    @test stats.num_constraints == 5
    @test stats.num_linear_constraints == 3
    @test stats.num_linear_equalities == 1
    @test stats.num_linear_inequalities == 2
    @test stats.num_quadratic_constraints == 2
    @test stats.num_quadratic_equalities == 1
    @test stats.num_quadratic_inequalities == 1
    @test stats.num_equality_constraints == 2
    @test stats.num_inequality_constraints == 3

    @test stats.num_linear_terms == 7
    @test stats.num_quadratic_terms == 4
    @test stats.num_diagonal_quadratic_terms == 2
    @test stats.num_bilinear_terms == 2
    @test stats.num_nonzeros == 11
    @test stats.num_distinct_coefficients == 10
    @test stats.num_constraints_with_fractional_coefficients == 1

    @test stats.min_abs_coefficient == 1.5
    @test stats.max_abs_coefficient == 8.0
    @test isapprox(stats.avg_abs_coefficient, 45.5 / 11.0; atol = 1.0e-12)
    @test stats.median_abs_coefficient == 4.0
    @test stats.coefficient_abs_range == 6.5
    @test stats.num_unit_coefficients == 0
    @test stats.num_even_coefficients == 6
    @test stats.num_odd_coefficients == 4
    @test stats.num_coefficients_divisible_by_4 == 3
    @test stats.num_coefficients_2mod4 == 3

    @test stats.num_even_linear_coefficients == 3.0
    @test stats.num_even_diagonal_coefficients == 1.0
    @test stats.num_even_bilinear_coefficients == 2.0
    @test stats.num_divisible_by_4_linear_coefficients == 2.0
    @test stats.num_divisible_by_4_diagonal_coefficients == 0.0
    @test stats.num_constraints_all_linear_coefficients_even == 2
    @test stats.num_constraints_all_diagonal_linear_coefficients_even == 1
    @test stats.num_constraints_with_odd_diagonal_or_linear_coefficient == 1

    @test stats.max_abs_linear_coefficient == 7.0
    @test stats.avg_abs_linear_coefficient == 4.0
    @test stats.max_abs_diagonal_coefficient == 6.0
    @test stats.avg_abs_diagonal_coefficient == 3.75
    @test stats.max_abs_bilinear_coefficient == 8.0
    @test stats.avg_abs_bilinear_coefficient == 5.0

    @test stats.min_constraint_bound_width == 8.0
    @test stats.max_constraint_bound_width == 10.0
    @test stats.avg_constraint_bound_width == 9.0
    @test stats.median_constraint_bound_width == 9.0
    @test stats.min_constraint_lower_bound == -2.0
    @test stats.max_constraint_lower_bound == 10.0
    @test stats.avg_constraint_lower_bound == 3.75
    @test stats.min_constraint_upper_bound == 6.0
    @test stats.max_constraint_upper_bound == 10.0
    @test stats.avg_constraint_upper_bound == 8.6
    @test stats.max_abs_constraint_bound == 10.0
    @test isapprox(stats.avg_abs_constraint_bound, 62.0 / 9.0; atol = 1.0e-12)
    @test stats.median_abs_constraint_bound == 8.0
end

@testset "ModelStats counts canonical bilinear terms once after permutation" begin
    qe = PC.QuadExpr(
        StatsQuadTerm[(2.0, 1, 1), (4.0, 2, 3), (5.0, 3, 2)],
        StatsLinTerm[],
    )
    PC.remove_var!(qe, 1)
    con = PC.Constraint(model_stats_next_con_id(), qe, 0.0, 12.0)
    model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(
            2 => PC.IntVar(0.0, 2.0),
            3 => PC.IntVar(0.0, 2.0),
        ),
        [con],
        model_stats_empty_objective(),
        :min,
    )

    stats = PC.ModelStats(model)

    @test stats.num_quadratic_terms == 1
    @test stats.num_diagonal_quadratic_terms == 0
    @test stats.num_bilinear_terms == 1
    @test stats.num_nonzeros == 1
    @test stats.max_abs_coefficient == 9.0
end

@testset "ModelStats returns zeros for empty collections" begin
    model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(),
        PC.Constraint[],
        model_stats_empty_objective(),
        :min,
    )

    stats = PC.ModelStats(model)

    @test stats.num_variables == 0
    @test stats.num_constraints == 0
    @test stats.num_nonzeros == 0
    @test stats.max_domain_size == 0
    @test stats.avg_domain_size == 0.0
    @test stats.median_domain_size == 0.0
    @test stats.min_abs_coefficient == 0.0
    @test stats.max_abs_coefficient == 0.0
    @test stats.avg_abs_coefficient == 0.0
    @test stats.median_abs_coefficient == 0.0
    @test stats.coefficient_abs_range == 0.0
    @test stats.min_constraint_bound_width == 0.0
    @test stats.max_constraint_bound_width == 0.0
    @test stats.avg_constraint_bound_width == 0.0
    @test stats.median_constraint_bound_width == 0.0
    @test stats.min_constraint_lower_bound == 0.0
    @test stats.max_constraint_lower_bound == 0.0
    @test stats.avg_constraint_lower_bound == 0.0
    @test stats.min_constraint_upper_bound == 0.0
    @test stats.max_constraint_upper_bound == 0.0
    @test stats.avg_constraint_upper_bound == 0.0
    @test stats.max_abs_constraint_bound == 0.0
    @test stats.avg_abs_constraint_bound == 0.0
    @test stats.median_abs_constraint_bound == 0.0
end
