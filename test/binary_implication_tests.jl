using Test
import Graphs
import QIPresolve.PresolvingCore as PC

function binary_test_row(id, coefficients, lo, hi; quadratic = Tuple{Float64, Int, Int}[])
    linear = Tuple{Float64, Int}[(a, i) for (i, a) in enumerate(coefficients) if !iszero(a)]
    return PC.Constraint(id, PC.QuadExpr(quadratic, linear), Float64(lo), Float64(hi))
end

function binary_test_model(n, cons)
    objective = PC.QuadExpr(Tuple{Float64, Int, Int}[], Tuple{Float64, Int}[])
    return PC.QPModel(Dict(i => PC.IntVar(0.0, 1.0) for i in 1:n), cons, objective, :min)
end

# Inspect the seeded graph before SCC condensation, independently of propagation.
function seed_accepts(manager, point)
    manager.infeasible && return false
    value(v) = begin
        lit = manager.pos_to_lit[v]
        Bool(point[lit.vid]) ⊻ lit.neg
    end
    return all(e -> !value(Graphs.src(e)) || value(Graphs.dst(e)), Graphs.edges(manager.graph))
end

@testset "exact binary implications for every two-variable relation" begin
    points = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for mask in 0:15
        allowed = [isodd(mask >> i) for i in 0:3]
        # Distinct activities keep both variables in the row, even for unit clauses.
        c, y, x, xy = [i - 1 + (allowed[i] ? 0 : 8) for i in 1:4]
        row = binary_test_row(1, [x - c, y - c], -Inf, 3 - c;
            quadratic = [(Float64(xy - x - y + c), 1, 2)])
        manager = PC.PropagationManager(Int[])
        model = binary_test_model(2, [row])
        PC._seed_binary_implications!(manager, model)
        @test [seed_accepts(manager, p) for p in points] == allowed
        PC.update!(manager)
        PC._seed_binary_implications!(manager, model)
        @test !manager.graph_dirty
    end

    # 2^53 + 1 rounds down in Float64. The (1, 1) assignment exceeds the bound.
    big = 2.0^53
    model = binary_test_model(2, [binary_test_row(1, [big, 1], -Inf, big)])
    manager = PC.PropagationManager(Int[])
    PC._seed_binary_implications!(manager, model)
    @test [seed_accepts(manager, p) for p in points] == [true, true, true, false]
    model.vars[1] = PC.IntVar(0, 2)
    unseeded = PC._seed_binary_implications!(PC.PropagationManager(Int[]), model)
    @test isempty(unseeded.pos_to_lit)
end

@testset "pair projections of longer signed binary inequalities" begin
    points = collect(Iterators.product(0:1, 0:1, 0:1))
    for a in Iterators.product((-2, 1), (-1, 2), (-2, 1)), bound in -4:5,
            (lo, hi) in ((-Inf, bound), (bound, Inf), (bound - 2, bound))
        model = binary_test_model(3, [binary_test_row(1, a, lo, hi)])
        manager = PC._seed_binary_implications!(PC.PropagationManager(Int[]), model)
        activity(p) = sum(a .* p)
        for p in points
            # Each side independently must admit an extension of each pair.
            expected = all(((i, j),) -> begin
                    completions = [activity(q) for q in points if q[i] == p[i] && q[j] == p[j]]
                    minimum(completions) <= hi && maximum(completions) >= lo
                end, ((1, 2), (1, 3), (2, 3)))
            @test seed_accepts(manager, p) == expected
        end
    end
end

@testset "inequality-only fixings, contradictions and phase re-scans" begin
    # u ≤ z and u + z ≥ 1 force z. The next phase sees x + y ≤ 1,
    # which combines with x ≤ y to force x = 0.
    rows = [binary_test_row(1, [1, 1, 1], -Inf, 2),
        binary_test_row(2, [1, -1], -Inf, 0),
        binary_test_row(3, [0, 0, -1, 1], -Inf, 0),
        binary_test_row(4, [0, 0, 1, 1], 1, Inf)]
    for strategy in (:mod2_basic, :mod4_basic, :full)
        model = binary_test_model(4, deepcopy(rows))
        manager = PC.PropagationManager(Int[])
        post = PC.ParityPostsolver(keys(model.vars))
        @test PC.parity_presolve_phase!(model, manager, post; parity_strategy = strategy)
        @test haskey(model.vars, 1) && !haskey(model.vars, 3)
        @test PC.parity_presolve_phase!(model, manager, post; parity_strategy = strategy)
        @test !model.infeasible && Set(keys(model.vars)) == Set([2, 4])
        for y in 0:1, u in 0:1
            @test PC.postsolve(post, Dict(2 => y, 4 => u)) == Dict(1 => 0, 2 => y, 3 => 1, 4 => u)
        end
    end
    # Four nonparallel inequalities each exclude a different corner.
    impossible = [binary_test_row(1, [3, -2], -Inf, 1), binary_test_row(2, [-2, 3], -Inf, 1),
        binary_test_row(3, [2, 1], -Inf, 2), binary_test_row(4, [1, 1], 1, Inf)]
    @test PC.presolve!(binary_test_model(2, impossible); enable_residue = false).model.infeasible
    disabled = PC.presolve!(binary_test_model(4, deepcopy(rows)); enable_parity = false, enable_residue = false)
    @test length(disabled.model.vars) == 4
end

@testset "graph-only representatives and implication bridges between XOR blocks" begin
    # p₁ = ¬p₂, p₁ = p₃. Only 2 and 3 have XOR columns.
    model = binary_test_model(3, [binary_test_row(1, [1, 1], 1, 1),
        binary_test_row(2, [1, 0, -1], 0, 0)])
    manager = PC.PropagationManager([1])
    PC._seed_binary_implications!(manager, model)
    xor_row = PC.XorConstraintBuilder()
    for vid in [2, 3, 4]
        PC.add_par!(xor_row, vid)
    end
    blocks = PC._build_parity_blocks!([xor_row], manager)
    @test !haskey(blocks.owner, 1)
    @test PC._reduce_parity_blocks!(blocks, manager, :full)
    @test PC.fixed_value(manager, 4) === true

    x, y = PC.XorConstraintBuilder(), PC.XorConstraintBuilder()
    PC.add_par!(x, 2)
    PC.add_par!(y, 3)
    graph = PC.PropagationManager(Int[])
    PC._seed_binary_implications!(graph, model)
    @test length(PC._build_parity_blocks!([x, y], graph).models) == 1
    @test PC._reduce_parity_blocks!(PC._build_parity_blocks!([x], graph), graph, :full)
    @test PC.fixed_value(graph, 1) === true
    @test PC.fixed_value(graph, 3) === true
end

@testset "odd domain shifts preserve retained implication polarity" begin
    manager = PC.PropagationManager([1, 2])
    PC.add_implication!(manager, PC.VarLit(1, false), PC.VarLit(2, false))
    model = binary_test_model(2, PC.Constraint[])
    model.vars[1] = PC.IntVar(3, 4)
    PC.normalize!(model; propagator = manager)
    # The shifted value zero represents the original odd value three.
    PC.fix_var!(manager, 1, false)
    PC.update!(manager)
    @test !manager.infeasible && PC.fixed_value(manager, 2) === true
end

@testset "residue normalization exposes new binary implications" begin
    rows = [binary_test_row(1, [1], -Inf, 1), binary_test_row(2, [1, -1], -Inf, 0),
        binary_test_row(3, [1, 1], 1, Inf)]
    model = binary_test_model(2, rows)
    model.vars[1] = PC.IntVar(0, 2)
    manager = PC.PropagationManager(Int[])
    PC._seed_binary_implications!(manager, model)
    @test isempty(manager.pos_to_lit)
    post = PC.ParityPostsolver(keys(model.vars))
    pass = PC._residue_presolve_pass!(model, Int[], nothing, post; treewidth_threshold = 2,
        propagator = manager)
    @test pass.changed && !pass.tightened_to_equality
    @test PC.parity_presolve_phase!(model, manager, post)
    @test !model.infeasible && !haskey(model.vars, 2)
    @test PC.postsolve(post, Dict(1 => 0)) == Dict(1 => 0, 2 => 1)
end
