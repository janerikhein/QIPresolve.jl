using Test
import QIPresolve.PresolvingCore as PC

const ParityQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const ParityLinTerm = Tuple{Float64, PC.VarId}
const PARITY_NEXT_CON_ID = Ref(0)

parity_next_con_id() = (PARITY_NEXT_CON_ID[] += 1)

function parity_empty_objective()
    return PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[])
end

function parity_edge_matrix(n::Int, edges::Vector{Tuple{Int, Int}})
    mat = falses(n, n)
    for (i, j) in edges
        mat[i, j] = true
        mat[j, i] = true
    end
    return mat
end

@testset "count_builder_occurrences! ignores toggled-off builder entries" begin
    counts = Dict{PC.VarId, Int}()
    builder = PC.XorConstraintBuilder()

    PC.add_par!(builder, 7)
    PC.add_par!(builder, 7)
    PC.add_par!(builder, 3)

    PC.add_conj!(builder, 2, 5)
    PC.add_conj!(builder, 5, 2)
    PC.add_conj!(builder, 3, 5)

    PC.count_builder_occurrences!(counts, builder)

    @test counts == Dict(3 => 2, 5 => 1)
    @test !haskey(counts, 2)
    @test !haskey(counts, 7)
end

@testset "build_parity_model orders active vars by occurrence count and preserves row order" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(-1.0, 2.0),
        2 => PC.IntVar(-1.0, 2.0),
        3 => PC.IntVar(0.0, 1.0),
        4 => PC.IntVar(0.0, 1.0),
    )

    con1 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 2)]),
        0.0,
        0.0,
    )
    con2 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 3)]),
        0.0,
        0.0,
    )
    con3 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 2)], ParityLinTerm[]),
        0.0,
        0.0,
    )
    con4 = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 2, 3)], ParityLinTerm[]),
        0.0,
        0.0,
    )

    model = PC.QPModel(vars, [con1, con2, con3, con4], parity_empty_objective(), :min)
    parity_model = PC.build_parity_model(model)

    @test parity_model.pos_to_var_id == [2, 1, 3]
    @test parity_model.var_id_to_pos == Dict(2 => 1, 1 => 2, 3 => 3)
    @test !haskey(parity_model.var_id_to_pos, 4)
    @test length(parity_model.cons) == 6

    @test parity_model.cons[1].meta.is_pure_xor
    @test parity_model.cons[1].par == BitVector([1, 0, 0])
    @test !parity_model.cons[1].rhs

    @test parity_model.cons[2].meta.is_pure_xor
    @test parity_model.cons[2].par == BitVector([0, 1, 1])
    @test !parity_model.cons[2].rhs

    @test parity_model.cons[3].meta.is_pure_xor
    @test parity_model.cons[3].par == falses(3)
    @test parity_model.cons[3].conj === nothing

    @test !parity_model.cons[4].meta.is_pure_xor
    @test parity_model.cons[4].par == falses(3)
    @test parity_model.cons[4].conj == parity_edge_matrix(3, [(1, 2)])

    @test parity_model.cons[5].meta.is_pure_xor
    @test parity_model.cons[5].par == falses(3)
    @test parity_model.cons[5].conj === nothing

    @test !parity_model.cons[6].meta.is_pure_xor
    @test parity_model.cons[6].par == falses(3)
    @test parity_model.cons[6].conj == parity_edge_matrix(3, [(1, 3)])
end


@testset "parity XOR branch expands conjunction terms using pivot substitutions" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
        3 => PC.IntVar(0.0, 1.0),
        4 => PC.IntVar(0.0, 1.0),
    )

    xor_con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[], ParityLinTerm[(1.0, 1), (1.0, 2), (1.0, 3)]),
        1.0,
        1.0,
    )
    mixed_con = PC.Constraint(
        parity_next_con_id(),
        PC.QuadExpr(ParityQuadTerm[(2.0, 1, 4)], ParityLinTerm[(2.0, 4)]),
        0.0,
        0.0,
    )

    model = PC.QPModel(vars, [xor_con, mixed_con], parity_empty_objective(), :min)
    parity_model = PC.build_parity_model(model)
    propagator = PC.PropagationManager(parity_model.pos_to_var_id)

    PC.reformulate_bipartite_cons!(parity_model)
    PC.propagate!(parity_model, propagator)
    PC.gauss_jordan_xor!(parity_model)
    PC.propagate!(parity_model, propagator)
    PC.substitute_pivots_in_conjunctive_terms!(parity_model)

    @test parity_model.pos_to_var_id == [1, 2, 3, 4]
    @test length(parity_model.cons) == 3
    @test parity_model.pivots[1] == (1, nothing)
    @test parity_model.cons[1].par == BitVector([1, 1, 1, 0])
    @test parity_model.cons[1].rhs
    @test parity_model.cons[2].par == falses(4)
    @test parity_model.cons[2].conj == parity_edge_matrix(4, [(2, 3)])
    @test !parity_model.cons[2].rhs
    @test parity_model.cons[3].par == falses(4)
    @test parity_model.cons[3].conj == parity_edge_matrix(4, [(2, 4), (3, 4)])
    @test !parity_model.cons[3].rhs
end
