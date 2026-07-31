using Test
using Random
using QIPresolve
import QIPresolve.PresolvingCore as PC

const DPQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const DPLinTerm = Tuple{Float64, PC.VarId}

function residue_bits(modulus::Int, residues::AbstractVector{Int})
    bits = falses(modulus)
    for residue in residues
        bits[mod(residue, modulus) + 1] = true
    end
    return bits
end

function conditioned_index(set::PC.ConditionedResidueSet, values::Dict{PC.VarId, Int})
    index = 1
    for (position, var_id) in enumerate(set.var_ids)
        index += values[var_id] * set.strides[position]
    end
    return index
end

function brute_component_residues(
        component::PC.NonSingleton,
        modulus::Int,
        var_ubs::Dict{PC.VarId, Int},
    )
    result = falses(modulus)
    values = Dict{PC.VarId, Int}()

    function visit(position::Int)
        if position > length(component.pos_to_var_id)
            total = 0
            for (var_id, coefficient) in component.lin_coeffs
                total += coefficient * values[var_id]
            end
            for ((first_id, second_id), coefficient) in component.quad_coeffs
                total += coefficient * values[first_id] * values[second_id]
            end
            result[mod(total, modulus) + 1] = true
            return nothing
        end

        var_id = component.pos_to_var_id[position]
        for value in 0:var_ubs[var_id]
            values[var_id] = value
            visit(position + 1)
        end
        return nothing
    end

    visit(1)
    return result
end

function capped_bounds(
        bounds::Dict{PC.VarId, PC.IntVar},
        modulus::Int,
        ids::AbstractVector{PC.VarId},
    )
    return Dict(var_id => min(trunc(Int, bounds[var_id].ub), modulus - 1) for var_id in ids)
end

function brute_linear_convolution(
        residues::BitVector,
        coeff::Int,
        ub::Int,
        modulus::Int,
    )
    output = falses(modulus)
    for residue in 0:(modulus - 1)
        residues[residue + 1] || continue
        for value in 0:min(ub, modulus - 1)
            output[mod(residue + coeff * value, modulus) + 1] = true
        end
    end
    return output
end

function brute_quad_singleton_residues(
        component::PC.QuadSingleton,
        ub::Int,
        modulus::Int,
    )
    output = falses(modulus)
    for value in 0:min(ub, modulus - 1)
        residue = component.quad_coeff * value * value + component.lin_coeff * value
        output[mod(residue, modulus) + 1] = true
    end
    return output
end

function brute_constraint_residues(
        con::PC.Constraint,
        modulus::Int,
        bounds::Dict{PC.VarId, PC.IntVar},
    )
    result = falses(modulus)
    var_ids = sort!(collect(PC.vars(con.qe)))
    values = Dict{PC.VarId, Int}()

    function visit(position::Int)
        if position > length(var_ids)
            total = 0
            for var_id in var_ids
                total += trunc(Int, PC.get_lin_coeff(con.qe, var_id)) * values[var_id]
            end
            for (i, first_id) in enumerate(var_ids)
                for second_id in @view var_ids[i:end]
                    coefficient = PC.get_quad_coeff(con.qe, first_id, second_id)
                    total += trunc(Int, coefficient) * values[first_id] * values[second_id]
                end
            end
            result[mod(total, modulus) + 1] = true
            return nothing
        end

        var_id = var_ids[position]
        ub = min(trunc(Int, bounds[var_id].ub), modulus - 1)
        for value in 0:ub
            values[var_id] = value
            visit(position + 1)
        end
        return nothing
    end

    visit(1)
    return result
end

@testset "ConditionedResidueSet construction and conditioning" begin
    set = PC.ConditionedResidueSet(5)

    @test set.modulus == 5
    @test isempty(set.var_ids)
    @test isempty(set.var_ubs)
    @test isempty(set.strides)
    @test set.residues == [residue_bits(5, [0])]
    @test set.local_contribution == [0]
    @test_throws ArgumentError PC.ConditionedResidueSet(0)

    PC.condition_on!(set, 10, 1)
    @test set.var_ids == [10]
    @test set.var_ubs == [1]
    @test set.strides == [1]
    @test set.residues == [residue_bits(5, [0]), residue_bits(5, [0])]

    set.residues[1][2] = true
    @test !set.residues[2][2]

    PC.condition_on!(set, 5, 2)
    @test set.var_ids == [5, 10]
    @test set.var_ubs == [2, 1]
    @test set.strides == [2, 1]
    @test length(set.residues) == 6
    for value in 0:2
        @test set.residues[1 + value * 2] == residue_bits(5, [0, 1])
        @test set.residues[2 + value * 2] == residue_bits(5, [0])
    end

    @test_throws ArgumentError PC.condition_on!(set, 10, 1)
    @test_throws ArgumentError PC.condition_on!(set, 20, 5)
end

@testset "local contribution is lex indexed after Gray-code enumeration" begin
    set = PC.ConditionedResidueSet(7)
    PC.condition_on!(set, 2, 1)
    PC.condition_on!(set, 3, 2)
    PC.condition_on!(set, 1, 2)

    quad_terms = Dict{Tuple{PC.VarId, PC.VarId}, Int}(
        (2, 2) => 4,
        (1, 2) => 5,
        (2, 3) => 6,
    )
    lin_terms = Dict{PC.VarId, Int}(2 => 3)

    PC.populate_local_contribution!(
        set,
        2,
        quad_terms,
        lin_terms,
        Set{PC.VarId}([3]),
    )

    @test set.var_ids == [1, 2, 3]
    @test set.strides == [6, 3, 1]
    for x1 in 0:2, x2 in 0:1, x3 in 0:2
        index = conditioned_index(set, Dict(1 => x1, 2 => x2, 3 => x3))
        expected = mod(4 * x2^2 + 3 * x2 + 5 * x1 * x2, 7)
        @test set.local_contribution[index] == expected
    end
end

@testset "fold_on shifts residues and merges over removed variable" begin
    set = PC.ConditionedResidueSet(5)
    PC.condition_on!(set, 1, 1)
    PC.condition_on!(set, 2, 1)

    set.residues = [
        residue_bits(5, [0]),
        residue_bits(5, [1]),
        residue_bits(5, [2]),
        residue_bits(5, [3]),
    ]
    set.local_contribution = [0, 1, 2, 3]

    PC.fold_on!(set, 2)

    @test set.var_ids == [1]
    @test set.var_ubs == [1]
    @test set.strides == [1]
    @test set.residues[1] == residue_bits(5, [0, 2])
    @test set.residues[2] == residue_bits(5, [1, 4])
    @test set.local_contribution == [0, 0]
end

@testset "join uses modular sumsets per conditioned assignment" begin
    left = PC.ConditionedResidueSet(5)
    right = PC.ConditionedResidueSet(5)
    PC.condition_on!(left, 1, 1)
    PC.condition_on!(right, 1, 1)

    left.residues = [residue_bits(5, [1, 3]), residue_bits(5, [2])]
    right.residues = [residue_bits(5, [0, 4]), residue_bits(5, [4])]

    PC.join!(left, right)

    @test left.residues[1] == residue_bits(5, [0, 1, 2, 3])
    @test left.residues[2] == residue_bits(5, [1])
    @test left.local_contribution == [0, 0]

    mismatch = PC.ConditionedResidueSet(5)
    PC.condition_on!(mismatch, 2, 1)
    @test_throws ArgumentError PC.join!(left, mismatch)
end

@testset "nonlinear component DP matches brute force" begin
    path_con = PC.Constraint(
        1,
        PC.QuadExpr(
            DPQuadTerm[
                (4.0, 1, 1),
                (2.0, 2, 2),
                (3.0, 1, 2),
                (5.0, 2, 3),
            ],
            DPLinTerm[(1.0, 1), (2.0, 2), (3.0, 3)],
        ),
        -100.0,
        100.0,
    )
    path_bounds = Dict(
        1 => PC.IntVar(0.0, 2.0),
        2 => PC.IntVar(0.0, 2.0),
        3 => PC.IntVar(0.0, 1.0),
    )

    path_results = PC.compute_nonlinear_residue_sets(7, path_con, path_bounds)
    @test length(path_results) == 1
    path_expected = brute_component_residues(
        path_results[1].component,
        7,
        capped_bounds(path_bounds, 7, path_results[1].component.pos_to_var_id),
    )
    @test path_results[1].residues == path_expected
    @test path_results[1].saturated == all(path_expected)

    star_con = PC.Constraint(
        2,
        PC.QuadExpr(
            DPQuadTerm[
                (2.0, 30, 10),
                (3.0, 30, 20),
                (5.0, 30, 40),
                (7.0, 30, 50),
                (4.0, 30, 30),
            ],
            DPLinTerm[
                (1.0, 10),
                (2.0, 20),
                (3.0, 30),
                (4.0, 40),
                (5.0, 50),
            ],
        ),
        -100.0,
        100.0,
    )
    star_bounds = Dict(
        10 => PC.IntVar(0.0, 1.0),
        20 => PC.IntVar(0.0, 1.0),
        30 => PC.IntVar(0.0, 1.0),
        40 => PC.IntVar(0.0, 1.0),
        50 => PC.IntVar(0.0, 1.0),
    )

    star_component = only(PC.decompose(PC.InteractionGraph(star_con, 17)))
    @test any(action isa PC.JoinAction for action in PC.action_order(star_component))

    star_results = PC.compute_nonlinear_residue_sets(17, star_con, star_bounds)
    @test length(star_results) == 1
    star_expected = brute_component_residues(
        star_results[1].component,
        17,
        capped_bounds(star_bounds, 17, star_results[1].component.pos_to_var_id),
    )
    @test star_results[1].residues == star_expected
    @test star_results[1].saturated == all(star_expected)
end

@testset "top-level DP validates bounds and caps locally" begin
    con = PC.Constraint(
        3,
        PC.QuadExpr(DPQuadTerm[(1.0, 1, 2)], DPLinTerm[]),
        -10.0,
        10.0,
    )

    @test_throws ArgumentError PC.compute_nonlinear_residue_sets(
        0,
        con,
        Dict(1 => PC.IntVar(0.0, 1.0), 2 => PC.IntVar(0.0, 1.0)),
    )
    @test_throws ArgumentError PC.compute_nonlinear_residue_sets(
        3,
        con,
        Dict(1 => PC.IntVar(0.0, 1.0)),
    )
    @test_throws ArgumentError PC.compute_nonlinear_residue_sets(
        3,
        con,
        Dict(1 => PC.IntVar(1.0, 2.0), 2 => PC.IntVar(0.0, 1.0)),
    )
    @test_throws ArgumentError PC.compute_nonlinear_residue_sets(
        3,
        con,
        Dict(1 => PC.IntVar(0.0, Inf), 2 => PC.IntVar(0.0, 1.0)),
    )
    @test_throws ArgumentError PC.compute_nonlinear_residue_sets(
        3,
        con,
        Dict(1 => PC.IntVar(0.0, 1.5), 2 => PC.IntVar(0.0, 1.0)),
    )

    bounds = Dict(
        1 => PC.IntVar(0.0, 9.0),
        2 => PC.IntVar(0.0, 9.0),
    )
    results = PC.compute_nonlinear_residue_sets(3, con, bounds)

    @test bounds[1].ub == 9.0
    @test bounds[2].ub == 9.0
    @test length(results) == 1
    expected = brute_component_residues(
        results[1].component,
        3,
        Dict(1 => 2, 2 => 2),
    )
    @test results[1].residues == expected
end

@testset "treewidth threshold saturates wide nonlinear components" begin
    path_con = PC.Constraint(
        4,
        PC.QuadExpr(
            DPQuadTerm[(3.0, 1, 2), (5.0, 2, 3)],
            DPLinTerm[(1.0, 1), (2.0, 2), (3.0, 3)],
        ),
        -100.0,
        100.0,
    )
    path_bounds = Dict(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
        3 => PC.IntVar(0.0, 1.0),
    )
    path_component = only(PC.decompose(PC.InteractionGraph(path_con, 7)))
    path_td = PC.minimum_degree_tree_decomposition(path_component)
    @test PC.tree_decomposition_width(path_td) == 1

    path_results = PC.compute_nonlinear_residue_sets(
        7,
        path_con,
        path_bounds;
        treewidth_threshold = 1,
    )
    path_expected = brute_component_residues(
        path_results[1].component,
        7,
        capped_bounds(path_bounds, 7, path_results[1].component.pos_to_var_id),
    )
    @test path_results[1].residues == path_expected
    @test !path_results[1].saturated

    triangle_con = PC.Constraint(
        5,
        PC.QuadExpr(
            DPQuadTerm[(2.0, 1, 2), (2.0, 1, 3), (2.0, 2, 3)],
            DPLinTerm[],
        ),
        1.0,
        7.0,
    )
    triangle_bounds = Dict(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
        3 => PC.IntVar(0.0, 1.0),
    )
    triangle_component = only(PC.decompose(PC.InteractionGraph(triangle_con, 4)))
    triangle_td = PC.minimum_degree_tree_decomposition(triangle_component)
    @test PC.tree_decomposition_width(triangle_td) == 2

    saturated_result = PC._compute_achievable_residues(
        4,
        triangle_con,
        triangle_bounds;
        treewidth_threshold = 1,
    )
    @test saturated_result.saturated
    @test saturated_result.residues == trues(4)

    exact_result = PC._compute_achievable_residues(
        4,
        triangle_con,
        triangle_bounds;
        treewidth_threshold = 2,
    )
    @test exact_result.residues == brute_constraint_residues(triangle_con, 4, triangle_bounds)
    @test !exact_result.saturated

    skipped_model = PC.QPModel(
        deepcopy(triangle_bounds),
        [deepcopy(triangle_con)],
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[]),
        :min,
    )
    PC.residue_presolve!(
        skipped_model,
        :divisor_free;
        threshold = 4,
        treewidth_threshold = 1,
    )
    @test !skipped_model.infeasible
    @test skipped_model.cons[1].lhs == triangle_con.lhs
    @test skipped_model.cons[1].rhs == triangle_con.rhs

    processed_model = PC.QPModel(
        deepcopy(triangle_bounds),
        [deepcopy(triangle_con)],
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[]),
        :min,
    )
    PC.residue_presolve!(
        processed_model,
        :divisor_free;
        threshold = 4,
        treewidth_threshold = 2,
    )
    @test !processed_model.infeasible
    @test processed_model.cons[1].lhs == 2.0
    @test processed_model.cons[1].rhs == 6.0
end

@testset "residue domain standardization shifts variables and objective" begin
    con = PC.Constraint(
        4,
        PC.QuadExpr(
            DPQuadTerm[(2.0, 1, 1), (4.0, 1, 2)],
            DPLinTerm[(3.0, 1), (-5.0, 2)];
            constant = 7.0,
        ),
        -10.0,
        50.0,
    )
    bounds = Dict(
        1 => PC.IntVar(-2.0, 3.0),
        2 => PC.IntVar(5.0, 7.0),
    )
    obj = PC.QuadExpr(
        DPQuadTerm[(1.0, 1, 2)],
        DPLinTerm[(2.0, 1), (3.0, 2)];
        constant = 11.0,
    )

    standardized = PC._standardize_residue_constraint(con, bounds, obj)
    lower_point = zeros(2)
    lower_point[1] = bounds[1].lb
    lower_point[2] = bounds[2].lb

    @test standardized.var_bounds == Dict(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(0.0, 2.0),
    )
    @test standardized.var_ubs == Dict(1 => 5, 2 => 2)
    @test standardized.constraint_shift == PC.eval_full(con.qe, lower_point)
    @test standardized.objective_shift == PC.eval_full(obj, lower_point) - obj.constant
    @test standardized.con.lhs == con.lhs - standardized.constraint_shift
    @test standardized.con.rhs == con.rhs - standardized.constraint_shift

    for x1 in -2:3, x2 in 5:7
        original = zeros(2)
        shifted = zeros(2)
        original[1] = x1
        original[2] = x2
        shifted[1] = x1 - bounds[1].lb
        shifted[2] = x2 - bounds[2].lb

        @test PC.eval_full(standardized.con.qe, shifted) + standardized.constraint_shift ==
            PC.eval_full(con.qe, original)
    end
end

@testset "residue modulus generation" begin
    @test PC._generate_residue_moduli(:small_primes, 1) == Int[]
    @test PC._generate_residue_moduli(:small_primes, 12) == [2, 3, 5, 7, 11]
    @test PC._generate_residue_moduli(:powers_of_two, 1) == Int[]
    @test PC._generate_residue_moduli(:powers_of_two, 20) == [2, 4, 8, 16]
    @test PC._generate_residue_moduli(:divisor_free, 1) == Int[]
    @test PC._generate_residue_moduli(:divisor_free, 2) == [2]
    @test PC._generate_residue_moduli(:divisor_free, 7) == [4, 5, 6, 7]
    @test PC._generate_residue_moduli(:divisor_free, 10) == [6, 7, 8, 9, 10]
    for threshold in 1:50
        moduli = PC._generate_residue_moduli(:divisor_free, threshold)
        @test issorted(moduli)
        @test all(2 <= modulus <= threshold for modulus in moduli)
        for (index, left) in enumerate(moduli)
            for right in @view moduli[(index + 1):end]
                @test left % right != 0
                @test right % left != 0
            end
        end
    end
    @test_throws ArgumentError PC._generate_residue_moduli(:unknown, 10)
end

@testset "linear singleton convolution matches brute force" begin
    explicit_cases = [
        (5, 2, 0, [0]),
        (8, 3, 4, [0, 5]),
        (11, -4, 10, [2, 7]),
        (13, 6, 18, [1, 3, 8]),
    ]

    for (modulus, coeff, ub, initial_residues) in explicit_cases
        actual = residue_bits(modulus, initial_residues)
        expected = brute_linear_convolution(copy(actual), coeff, ub, modulus)
        PC._convolve_linear_singleton!(actual, PC.LinSingleton(1, coeff), ub, modulus)
        @test actual == expected
    end

    rng = MersenneTwister(20260731)
    for _ in 1:100
        modulus = rand(rng, 2:17)
        coeff = rand(rng, -30:30)
        ub = rand(rng, 0:(2 * modulus + 3))
        initial = falses(modulus)
        for residue in 0:(modulus - 1)
            rand(rng, Bool) && (initial[residue + 1] = true)
        end
        initial[rand(rng, 1:modulus)] = true

        actual = copy(initial)
        PC._convolve_linear_singleton!(actual, PC.LinSingleton(1, coeff), ub, modulus)
        @test actual == brute_linear_convolution(initial, coeff, ub, modulus)
    end
end

@testset "quadratic singleton residues enumerate capped domain" begin
    for modulus in 2:13, lin_coeff in -3:3, quad_coeff in -4:4, ub in 0:(modulus + 2)
        quad_coeff == 0 && continue
        component = PC.QuadSingleton(1, lin_coeff, quad_coeff)
        @test PC._quad_singleton_residues(component, ub, modulus) ==
            brute_quad_singleton_residues(component, ub, modulus)
    end
end

@testset "global achievable residues merge all component classes" begin
    con = PC.Constraint(
        5,
        PC.QuadExpr(
            DPQuadTerm[(3.0, 2, 2), (4.0, 3, 4)],
            DPLinTerm[(2.0, 1), (1.0, 2), (5.0, 3), (-2.0, 4)],
        ),
        -100.0,
        100.0,
    )
    bounds = Dict(
        1 => PC.IntVar(0.0, 2.0),
        2 => PC.IntVar(0.0, 3.0),
        3 => PC.IntVar(0.0, 1.0),
        4 => PC.IntVar(0.0, 2.0),
    )
    modulus = 7

    result = PC._compute_achievable_residues(modulus, con, bounds)
    expected = brute_constraint_residues(con, modulus, bounds)
    @test result.residues == expected
    @test result.saturated == all(expected)
end

@testset "residue_presolve! tightens bounds in original coordinates" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(1.0, 4.0))
    con = PC.Constraint(
        6,
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[(2.0, 1)]),
        1.0,
        9.0,
    )
    obj = PC.QuadExpr(DPQuadTerm[(1.0, 1, 1)], DPLinTerm[(3.0, 1)]; constant = 4.0)
    model = PC.QPModel(vars, [con], deepcopy(obj), :min)

    @test QIPresolve.residue_presolve!(model, :powers_of_two; threshold = 2) === model
    @test !model.infeasible
    @test model.vars[1] == PC.IntVar(1.0, 4.0)
    @test model.cons[1].lhs == 2.0
    @test model.cons[1].rhs == 8.0
    @test model.obj_expr.constant == obj.constant
    @test PC.get_quad_coeff(model.obj_expr, 1, 1) == PC.get_quad_coeff(obj, 1, 1)
    @test PC.get_lin_coeff(model.obj_expr, 1) == PC.get_lin_coeff(obj, 1)
end

@testset "residue_presolve! reapplies cached moduli to fixed point" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 10.0))
    con = PC.Constraint(
        7,
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[(6.0, 1)]),
        1.0,
        20.0,
    )
    model = PC.QPModel(vars, [con], PC.QuadExpr(DPQuadTerm[], DPLinTerm[]), :min)

    PC.residue_presolve!(model, :small_primes; threshold = 3)
    @test !model.infeasible
    @test model.cons[1].lhs == 6.0
    @test model.cons[1].rhs == 18.0
end

@testset "residue_presolve! marks infeasible and skips unsupported domains" begin
    infeasible_model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 10.0)),
        [
            PC.Constraint(
                8,
                PC.QuadExpr(DPQuadTerm[], DPLinTerm[(2.0, 1)]),
                1.0,
                1.0,
            ),
        ],
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[]),
        :min,
    )
    PC.residue_presolve!(infeasible_model, :powers_of_two; threshold = 2)
    @test infeasible_model.infeasible

    nonfinite_model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, Inf)),
        [
            PC.Constraint(
                9,
                PC.QuadExpr(DPQuadTerm[], DPLinTerm[(2.0, 1)]),
                1.0,
                9.0,
            ),
        ],
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[]),
        :min,
    )
    PC.residue_presolve!(nonfinite_model, :powers_of_two; threshold = 2)
    @test !nonfinite_model.infeasible
    @test nonfinite_model.cons[1].lhs == 1.0
    @test nonfinite_model.cons[1].rhs == 9.0

    fractional_model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.5, 3.0)),
        [
            PC.Constraint(
                10,
                PC.QuadExpr(DPQuadTerm[], DPLinTerm[(2.0, 1)]),
                1.0,
                9.0,
            ),
        ],
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[]),
        :min,
    )
    PC.residue_presolve!(fractional_model, :powers_of_two; threshold = 2)
    @test !fractional_model.infeasible
    @test fractional_model.cons[1].lhs == 1.0
    @test fractional_model.cons[1].rhs == 9.0

    valid_model = PC.QPModel(
        Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0)),
        PC.Constraint[],
        PC.QuadExpr(DPQuadTerm[], DPLinTerm[]),
        :min,
    )
    @test_throws ArgumentError PC.residue_presolve!(valid_model, :unknown; threshold = 3)
end
