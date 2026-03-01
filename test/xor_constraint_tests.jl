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

@testset "XorConstraint constructors" begin
    n = 4
    pos_to_var, var_to_pos = identity_maps(n)
    par = BitVector([1, 0, 1, 0])
    conj = symmetric_bitmatrix(n, [(1, 2), (1, 3), (2, 4), (3, 4)])

    con = PC.XorConstraint(pos_to_var, var_to_pos, copy(par), copy(conj), true)
    @test con.pos_to_var == pos_to_var
    @test con.var_to_pos == var_to_pos
    @test con.par == par
    @test con.conj == conj
    @test con.rhs
    @test con.conj_deg == [2, 2, 2, 2]
    @test con.nnz == 6
    @test !con.is_pure_xor
    @test !con.has_changed

    pure = PC.XorConstraint(pos_to_var, var_to_pos, copy(par), false)
    @test pure.conj === nothing
    @test pure.conj_deg === nothing
    @test pure.nnz == 2
    @test pure.is_pure_xor
    @test !pure.has_changed
    @test !pure.rhs
end

@testset "update! recomputes cache and detects pure XOR" begin
    n = 3
    pos_to_var, var_to_pos = identity_maps(n)
    par = BitVector([1, 1, 0])
    conj = symmetric_bitmatrix(n, [(1, 2)])
    con = PC.XorConstraint(pos_to_var, var_to_pos, copy(par), copy(conj), false)

    con.conj .= false
    con.has_changed = true
    res = PC.update!(con)

    @test res == false
    @test con.has_changed == false
    @test con.is_pure_xor
    @test con.conj === nothing
    @test con.conj_deg === nothing
    @test con.nnz == 2
end

@testset "xor_con! toggles linear/conjunctive terms and rhs" begin
    n = 3
    pos_to_var, var_to_pos = identity_maps(n)
    conj1 = symmetric_bitmatrix(n, [(1, 2), (2, 3)])
    conj2 = symmetric_bitmatrix(n, [(2, 3)])
    con1 = PC.XorConstraint(pos_to_var, var_to_pos, BitVector([1, 0, 1]), conj1, false)
    con2 = PC.XorConstraint(pos_to_var, var_to_pos, BitVector([1, 1, 0]), conj2, true)

    out = PC.xor_con!(con1, con2)
    @test out === con1
    @test con1.par == BitVector([0, 1, 1])
    @test con1.conj == symmetric_bitmatrix(n, [(1, 2)])
    @test con1.rhs
    @test !con1.has_changed
    @test con1.nnz == 3
    @test con1.conj_deg == [1, 1, 0]

    con3 = PC.XorConstraint(pos_to_var, var_to_pos, BitVector([1, 0, 0]), false)
    con4 = PC.XorConstraint(pos_to_var, var_to_pos, BitVector([0, 1, 0]), true)
    PC.xor_con!(con3, con4; update = false)
    @test con3.par == BitVector([1, 1, 0])
    @test con3.rhs
    @test con3.has_changed
end

@testset "split_bipartite returns two XOR constraints for full bipartite conjs" begin
    n = 4
    pos_to_var, var_to_pos = identity_maps(n)
    par = falses(n)
    conj = symmetric_bitmatrix(n, [(1, 3), (1, 4), (2, 3), (2, 4)])
    con = PC.XorConstraint(pos_to_var, var_to_pos, par, conj, true)

    split = PC.split_bipartite(con)
    @test split !== nothing
    con1, con2 = split
    @test con1.is_pure_xor
    @test con2.is_pure_xor
    @test con1.rhs
    @test con2.rhs

    expected = Set([BitVector([1, 1, 0, 0]), BitVector([0, 0, 1, 1])])
    got = Set([con1.par, con2.par])
    @test got == expected

    not_split_rhs = PC.XorConstraint(pos_to_var, var_to_pos, par, conj, false)
    @test PC.split_bipartite(not_split_rhs) === nothing

    triangle = symmetric_bitmatrix(3, [(1, 2), (2, 3), (1, 3)])
    pos_to_var3, var_to_pos3 = identity_maps(3)
    non_bip = PC.XorConstraint(pos_to_var3, var_to_pos3, falses(3), triangle, true)
    @test PC.split_bipartite(non_bip) === nothing
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
    @test !con.is_pure_xor
    @test con.par == BitVector([0, 1, 0])
    @test con.conj == symmetric_bitmatrix(n, [(2, 3)])
    @test !con.rhs

    pure_builder = PC.XorConstraintBuilder()
    PC.add_par!(pure_builder, 3)
    pure = PC.build(pure_builder, n, pos_to_var, var_to_pos)
    @test pure.is_pure_xor
    @test pure.par == BitVector([0, 0, 1])
    @test pure.conj === nothing
end
