using Test
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
