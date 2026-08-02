using Test
import QIPresolve
import QIPresolve.PresolvingCore as PC

const PresolveQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const PresolveLinTerm = Tuple{Float64, PC.VarId}
const PRESOLVE_NEXT_CON_ID = Ref(0)

presolve_next_con_id() = (PRESOLVE_NEXT_CON_ID[] += 1)
presolve_empty_objective() = PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[])

@testset "presolve! returns a result with usable postsolve data" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0))
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[(1.0, 1)]),
        1.0,
        1.0,
    )
    model = PC.QPModel(vars, [con], presolve_empty_objective(), :min)

    result = QIPresolve.presolve!(
        model;
        residue_strategy = :powers_of_two,
        residue_threshold = 2,
    )

    @test result isa QIPresolve.PresolveResult
    @test result.model === model
    @test isempty(model.vars)
    @test QIPresolve.postsolve(result.postsolver, Dict{PC.VarId, Float64}()) ==
        Dict{PC.VarId, Float64}(1 => 1.0)
end

@testset "presolve! runs the first residue pass without parity changes" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
    )
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[(2.0, 1), (3.0, 2)]),
        1.0,
        4.0,
    )
    model = PC.QPModel(vars, [con], presolve_empty_objective(), :min)

    result = QIPresolve.presolve!(model)

    @test result.model === model
    @test !model.infeasible
    @test model.cons[1].lhs == 2.0
    @test model.cons[1].rhs == 3.0
    @test model.vars == vars
end

@testset "presolve! uses residue-created equalities to trigger parity" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 2.0))
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[(2.0, 1, 1)], PresolveLinTerm[]),
        1.0,
        2.0,
    )
    model = PC.QPModel(vars, [con], presolve_empty_objective(), :min)

    result = QIPresolve.presolve!(
        model;
        residue_strategy = :powers_of_two,
        residue_threshold = 2,
    )

    @test !model.infeasible
    @test isempty(model.vars)
    @test QIPresolve.postsolve(result.postsolver, Dict{PC.VarId, Float64}()) ==
        Dict{PC.VarId, Float64}(1 => 1.0)
end

@testset "residue_presolve_constraint! reports non-equality tightening" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(1.0, 4.0))
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[(2.0, 1)]),
        1.0,
        9.0,
    )
    model = PC.QPModel(vars, [con], presolve_empty_objective(), :min)

    stats = PC.residue_presolve_constraint!(model, con, [2])

    @test stats.changed
    @test !stats.tightened_to_equality
    @test stats.processed
    @test !stats.infeasible
    @test con.lhs == 2.0
    @test con.rhs == 8.0
end

@testset "filtered residue passes only touch candidate constraints" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(1.0, 4.0),
        2 => PC.IntVar(1.0, 4.0),
    )
    con1 = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[(2.0, 1)]),
        1.0,
        9.0,
    )
    con2 = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[(2.0, 2)]),
        1.0,
        9.0,
    )
    model = PC.QPModel(vars, [con1, con2], presolve_empty_objective(), :min)

    stats = PC._residue_presolve_pass!(
        model,
        [2],
        Set([con1.id]);
        treewidth_threshold = 2,
    )

    @test stats.changed
    @test stats.processed == 1
    @test stats.processed_constraint_ids == [con1.id]
    @test con1.lhs == 2.0
    @test con1.rhs == 8.0
    @test con2.lhs == 1.0
    @test con2.rhs == 9.0
end

@testset "parity snapshots mark newly added constraints as residue candidates" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 5.0),
        2 => PC.IntVar(1.0, 6.0),
    )
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[(1.0, 1), (1.0, 2)]),
        4.0,
        4.0,
    )
    model = PC.QPModel(vars, [con], presolve_empty_objective(), :min)
    propagator = PC.PropagationManager(PC.VarId[])
    postsolver = PC.ParityPostsolver(keys(vars))

    stats = PC._parity_presolve_fixed_point!(model, propagator, postsolver)
    new_con_ids = setdiff(Set(con.id for con in model.cons), Set([con.id]))

    @test stats.domains_changed
    @test !isempty(new_con_ids)
    @test all(con_id in stats.coefficient_changed_constraint_ids for con_id in new_con_ids)
end
