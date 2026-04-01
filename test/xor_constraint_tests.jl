using Test
import QIPresolve.PresolvingCore as PC

function identity_maps(n::Int)
    pos_to_var = collect(1:n)
    var_to_pos = Dict(i => i for i in 1:n)
    return pos_to_var, var_to_pos
end

function symmetric_bitmatrix(n::Int, edges::Vector{Tuple{Int, Int}})
    mat = falses(n, n)
    for (i, j) in edges
        mat[i, j] = true
        mat[j, i] = true
    end
    return mat
end

function constraint_holds(con::PC.XorConstraint, x::AbstractVector{Bool})
    lhs = false
    n = length(con.par)

    for i in 1:n
        con.par[i] || continue
        lhs = xor(lhs, x[i])
    end

    if con.conj !== nothing
        for i in 1:n
            xi = x[i]
            for j in (i + 1):n
                con.conj[i, j] || continue
                lhs = xor(lhs, xi & x[j])
            end
        end
    end

    return lhs == con.rhs
end

function foreach_assignment_except(f::Function, n::Int, skip_idx::Int)
    for bits in 0:((1 << (n - 1)) - 1)
        x = falses(n)
        free_idx = 0
        for i in 1:n
            i == skip_idx && continue
            x[i] = ((bits >> free_idx) & 1) == 1
            free_idx += 1
        end
        f(x)
    end
    return nothing
end

function assert_single_substitution_equivalent(
    original::PC.XorConstraint,
    transformed::PC.XorConstraint,
    var_idx::Int,
    subst_idx::Int,
    neg::Bool,
)
    n = length(original.par)
    foreach_assignment_except(n, var_idx) do x
        x[var_idx] = xor(x[subst_idx], neg)
        @test constraint_holds(transformed, x) == constraint_holds(original, x)
    end
end

function assert_mask_substitution_equivalent(
    original::PC.XorConstraint,
    transformed::PC.XorConstraint,
    var_idx::Int,
    subst_mask::BitVector,
    neg::Bool,
)
    n = length(original.par)
    foreach_assignment_except(n, var_idx) do x
        subst_val = neg
        for j in eachindex(subst_mask)
            subst_mask[j] || continue
            subst_val = xor(subst_val, x[j])
        end
        x[var_idx] = subst_val
        @test constraint_holds(transformed, x) == constraint_holds(original, x)
    end
end

@testset "XorConstraint constructors" begin
    par = BitVector([1, 0, 1, 0])
    conj = symmetric_bitmatrix(4, [(1, 2), (1, 3), (2, 4), (3, 4)])

    con = PC.XorConstraint(copy(par), copy(conj), true)
    @test con.par == par
    @test con.conj == conj
    @test con.rhs
    @test con.meta.nnz_par == 2
    @test con.meta.nnz_conj == 4
    @test !con.meta.is_pure_xor
    @test !con.meta.requires_update
    @test con.meta.requires_prop

    pure = PC.XorConstraint(copy(par), false)
    @test pure.conj === nothing
    @test pure.meta.nnz_par == 2
    @test pure.meta.nnz_conj == 0
    @test pure.meta.is_pure_xor
    @test !pure.meta.requires_update
    @test pure.meta.requires_prop
    @test !pure.rhs
end

@testset "update! recomputes cache and detects pure XOR" begin
    par = BitVector([1, 1, 0])
    conj = symmetric_bitmatrix(3, [(1, 2)])
    con = PC.XorConstraint(copy(par), copy(conj), false)

    con.conj .= false
    con.meta.requires_update = true
    con.meta.requires_prop = false
    res = PC.update!(con)

    @test res === con
    @test !con.meta.requires_update
    @test !con.meta.requires_prop
    @test con.meta.is_pure_xor
    @test con.conj === nothing
    @test con.meta.nnz_par == 2
    @test con.meta.nnz_conj == 0
end

@testset "xor_con! toggles linear/conjunctive terms and rhs" begin
    n = 3
    conj1 = symmetric_bitmatrix(n, [(1, 2), (2, 3)])
    conj2 = symmetric_bitmatrix(n, [(2, 3)])
    con1 = PC.XorConstraint(BitVector([1, 0, 1]), conj1, false)
    con2 = PC.XorConstraint(BitVector([1, 1, 0]), conj2, true)

    out = PC.xor_con!(con1, con2)
    @test out === con1
    @test con1.par == BitVector([0, 1, 1])
    @test con1.conj == symmetric_bitmatrix(n, [(1, 2)])
    @test con1.rhs
    @test con1.meta.requires_update
    @test con1.meta.requires_prop

    PC.update!(con1)
    @test con1.meta.nnz_par == 2
    @test con1.meta.nnz_conj == 1
    @test !con1.meta.is_pure_xor

    con3 = PC.XorConstraint(BitVector([1, 0, 0]), false)
    con4 = PC.XorConstraint(BitVector([0, 1, 0]), true)
    PC.xor_con!(con3, con4)
    @test con3.par == BitVector([1, 1, 0])
    @test con3.rhs
    @test con3.meta.requires_update
    @test con3.meta.requires_prop
end

@testset "split_bipartite returns two XOR constraints for full bipartite conjs" begin
    par = falses(4)
    conj = symmetric_bitmatrix(4, [(1, 3), (1, 4), (2, 3), (2, 4)])
    con = PC.XorConstraint(par, conj, true)

    split = PC.split_bipartite(con)
    @test split !== nothing
    con1, con2 = split
    @test con1.meta.is_pure_xor
    @test con2.meta.is_pure_xor
    @test con1.rhs
    @test con2.rhs

    expected = Set([BitVector([1, 1, 0, 0]), BitVector([0, 0, 1, 1])])
    got = Set([con1.par, con2.par])
    @test got == expected

    not_split_rhs = PC.XorConstraint(par, conj, false)
    @test PC.split_bipartite(not_split_rhs) === nothing

    triangle = symmetric_bitmatrix(3, [(1, 2), (2, 3), (1, 3)])
    non_bip = PC.XorConstraint(falses(3), triangle, true)
    @test PC.split_bipartite(non_bip) === nothing
end

@testset "substitute_var! with single target preserves semantics and invariants" begin
    original = PC.XorConstraint(
        BitVector([1, 0, 0, 1]),
        symmetric_bitmatrix(4, [(1, 2), (1, 3), (2, 4)]),
        false,
    )
    con = deepcopy(original)

    out = PC.substitute_var!(con, 1, 2, false)

    @test out === con
    @test con.par == BitVector([0, 0, 0, 1])
    @test con.conj == symmetric_bitmatrix(4, [(2, 3), (2, 4)])
    @test !con.rhs
    @test all(.!con.conj[:, 1])
    @test all(.!con.conj[1, :])
    @test con.conj == permutedims(con.conj)
    @test con.meta.requires_update
    @test con.meta.requires_prop

    assert_single_substitution_equivalent(original, con, 1, 2, false)

    original_neg = PC.XorConstraint(
        BitVector([1, 1, 0, 0]),
        symmetric_bitmatrix(4, [(1, 2), (1, 3), (2, 4)]),
        false,
    )
    con_neg = deepcopy(original_neg)

    PC.substitute_var!(con_neg, 1, 2, true)

    @test con_neg.par == BitVector([0, 0, 1, 0])
    @test con_neg.conj == symmetric_bitmatrix(4, [(2, 3), (2, 4)])
    @test con_neg.rhs
    @test all(.!con_neg.conj[:, 1])
    @test all(.!con_neg.conj[1, :])
    @test con_neg.conj == permutedims(con_neg.conj)
    @test con_neg.meta.requires_update
    @test con_neg.meta.requires_prop

    assert_single_substitution_equivalent(original_neg, con_neg, 1, 2, true)
end

@testset "substitute_var! with substitution masks preserves semantics and invariants" begin
    original = PC.XorConstraint(
        BitVector([1, 0, 1, 0, 0]),
        symmetric_bitmatrix(5, [(1, 2), (1, 4), (1, 5)]),
        false,
    )
    con = deepcopy(original)
    subst_mask = BitVector([0, 1, 0, 1, 0])

    out = PC.substitute_var!(con, 1, subst_mask, true)

    @test out === con
    @test con.par == BitVector([0, 1, 1, 1, 1])
    @test con.conj == symmetric_bitmatrix(5, [(2, 5), (4, 5)])
    @test con.rhs
    @test all(.!con.conj[:, 1])
    @test all(.!con.conj[1, :])
    @test con.conj == permutedims(con.conj)
    @test con.meta.requires_update
    @test con.meta.requires_prop

    assert_mask_substitution_equivalent(original, con, 1, subst_mask, true)

    original_const = PC.XorConstraint(
        BitVector([1, 0, 1, 0]),
        symmetric_bitmatrix(4, [(1, 2), (1, 4), (2, 3)]),
        true,
    )
    con_const = deepcopy(original_const)
    empty_mask = falses(4)

    PC.substitute_var!(con_const, 1, empty_mask, true)

    @test con_const.par == BitVector([0, 1, 1, 1])
    @test con_const.conj == symmetric_bitmatrix(4, [(2, 3)])
    @test !con_const.rhs
    @test all(.!con_const.conj[:, 1])
    @test all(.!con_const.conj[1, :])
    @test con_const.conj == permutedims(con_const.conj)
    @test con_const.meta.requires_update
    @test con_const.meta.requires_prop

    assert_mask_substitution_equivalent(original_const, con_const, 1, empty_mask, true)
end

@testset "XorConstraintBuilder toggling and build" begin
    n = 3
    pos_to_var, var_to_pos = identity_maps(n)

    builder = PC.XorConstraintBuilder()
    @test PC.add_par!(builder, 1)
    @test !PC.add_par!(builder, 1)
    @test PC.add_par!(builder, 2)

    @test PC.add_conj!(builder, 3, 1)
    @test !PC.add_conj!(builder, 1, 3)
    @test PC.add_conj!(builder, 2, 3)
    @test haskey(builder.conj, (1, 3))
    @test haskey(builder.conj, (2, 3))

    @test PC.negate!(builder)
    @test !PC.negate!(builder)

    con = PC.build(builder, n, pos_to_var, var_to_pos)
    @test !con.meta.is_pure_xor
    @test con.par == BitVector([0, 1, 0])
    @test con.conj == symmetric_bitmatrix(n, [(2, 3)])
    @test !con.rhs

    pure_builder = PC.XorConstraintBuilder()
    PC.add_par!(pure_builder, 3)
    pure = PC.build(pure_builder, n, pos_to_var, var_to_pos)
    @test pure.meta.is_pure_xor
    @test pure.par == BitVector([0, 0, 1])
    @test pure.conj === nothing
end
