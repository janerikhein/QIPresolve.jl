using Test
using Graphs: has_edge, nv, ne
import QIPresolve.PresolvingCore as PC

@testset "ParityPropagator constructor and mappings" begin
    prop = PC.ParityPropagator(3)

    @test prop.nvars == 3
    @test nv(prop.scc_graph) == 6
    @test ne(prop.scc_graph) == 0
    @test prop.requires_update == false
    @test length(prop.sccs) == 6
    @test PC.ImpGraph === PC.ParityPropagator

    for vid in 1:3
        pos = PC.lit_to_pos(prop, vid, false)
        npos = PC.lit_to_pos(prop, vid, true)
        @test pos == vid
        @test npos == vid + 3
        @test PC.pos_to_lit(prop, pos) == (vid, false)
        @test PC.pos_to_lit(prop, npos) == (vid, true)

        @test count(prop.sccs[pos].lits) == 1
        @test prop.sccs[pos].lits[vid]
        @test count(prop.sccs[pos].lits_neg) == 0

        @test count(prop.sccs[npos].lits_neg) == 1
        @test prop.sccs[npos].lits_neg[vid]
        @test count(prop.sccs[npos].lits) == 0
    end

    @test_throws ArgumentError PC.ParityPropagator(-1)
end

@testset "ParityModel initializes propagator with model variable order" begin
    var_to_pos = Dict(10 => 1, 7 => 2, 42 => 3)
    pos_to_var = [10, 7, 42]
    model = PC.ParityModel(var_to_pos, pos_to_var, PC.XorConstraint[])

    @test model.propagator.var_id_to_var_pos == var_to_pos
    @test model.propagator.var_pos_to_var_id == pos_to_var
    @test PC.lit_to_pos(model.propagator, 42, false) == 3
    @test PC.pos_to_lit(model.propagator, 2) == (7, false)
end

@testset "Union-find helpers" begin
    prop = PC.ParityPropagator(3)

    @test PC.union_scc!(prop, 2, 3) == 2
    @test PC.get_scc_repr_pos!(prop, 3) == 2

    @test PC.union_scc!(prop, 1, 3) == 1
    @test PC.get_scc_repr_pos!(prop, 2) == 1
    @test PC.get_scc_repr_pos!(prop, 3) == 1
    @test prop.parent_pos[2] == 1
    @test prop.parent_pos[3] == 1
end

@testset "add_edge!, fix_var!, and update!" begin
    prop = PC.ParityPropagator(2)

    @test PC.add_edge!(prop, 1, false, 2, false)
    @test prop.requires_update
    @test !PC.add_edge!(prop, 1, false, 2, false)
    @test prop.requires_update

    @test PC.fix_var!(prop, 1, true)
    @test has_edge(prop.scc_graph, 3, 1) # ¬x1 -> x1
    @test PC.fix_var!(prop, 2, false)
    @test has_edge(prop.scc_graph, 2, 4) # x2 -> ¬x2

    @test PC.update!(prop) === prop
    @test prop.requires_update == false
    @test !PC.fix_var!(prop, 1, true) # already present
end

@testset "update! condenses SCCs and mappings" begin
    prop = PC.ParityPropagator(2)
    @test PC.update!(prop) === prop # no-op

    PC.add_edge!(prop, 1, false, 2, false)
    PC.add_edge!(prop, 2, false, 1, false)
    @test prop.requires_update

    PC.update!(prop)
    @test !prop.requires_update

    repr1 = PC.get_scc_repr_pos!(prop, 1)
    repr2 = PC.get_scc_repr_pos!(prop, 2)
    @test repr1 == repr2
    @test nv(prop.scc_graph) == 3

    node_idx = prop.var_pos_to_node_idx[repr1]
    merged_scc = prop.sccs[node_idx]
    @test merged_scc.lits[1]
    @test merged_scc.lits[2]
    @test !any(merged_scc.lits_neg)
end

@testset "propagate! basic and forced propagation" begin
    # No implications, no forced literals.
    prop0 = PC.ParityPropagator(2)
    feasible0, t0, f0 = PC.propagate!(prop0)
    @test feasible0
    @test !any(t0)
    @test !any(f0)

    # Force x1 = true, then propagate x1 -> x2 -> ¬x3.
    prop = PC.ParityPropagator(3)
    PC.fix_var!(prop, 1, true)
    PC.add_edge!(prop, 1, false, 2, false)
    PC.add_edge!(prop, 2, false, 3, true)

    feasible, lit_true, lit_false = PC.propagate!(prop)
    @test feasible

    @test lit_true[1]   # x1
    @test lit_true[2]   # x2
    @test lit_true[6]   # ¬x3

    @test lit_false[4]  # ¬x1
    @test lit_false[5]  # ¬x2
    @test lit_false[3]  # x3
end

@testset "propagate! infeasibility detection" begin
    # Literal and its negation in same SCC.
    prop1 = PC.ParityPropagator(1)
    PC.add_edge!(prop1, 1, false, 1, true)
    PC.add_edge!(prop1, 1, true, 1, false)
    feasible1, _, _ = PC.propagate!(prop1)
    @test !feasible1

    # True literal implies a literal that is forced false.
    prop2 = PC.ParityPropagator(2)
    PC.fix_var!(prop2, 1, true)
    PC.fix_var!(prop2, 2, true)
    PC.add_edge!(prop2, 1, false, 2, true) # x1 -> ¬x2 while x2 is true
    feasible2, _, _ = PC.propagate!(prop2)
    @test !feasible2
end
