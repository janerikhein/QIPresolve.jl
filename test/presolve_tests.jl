using Test
using JuMP: backend
import QIPresolve
import QIPresolve.PresolvingCore as PC
using QIPresolve.InstanceGeneration: generate_random_qip_model

const PresolveQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const PresolveLinTerm = Tuple{Float64, PC.VarId}
const PRESOLVE_NEXT_CON_ID = Ref(0)

presolve_next_con_id() = (PRESOLVE_NEXT_CON_ID[] += 1)
presolve_empty_objective() = PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[])

@testset "presolve config exposes current default values" begin
    @test QIPresolve.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_STRATEGY == :divisor_free
    @test QIPresolve.PresolveConfig.DEFAULT_PRESOLVE_RESIDUE_THRESHOLD == 64
    @test QIPresolve.PresolveConfig.DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD == 2
    @test QIPresolve.PresolveConfig.DEFAULT_PRESOLVE_PARITY_STRATEGY == :full
end

@testset "presolve symmetrization doubles canonical constraint coefficients" begin
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[(3.0, 1, 2)], PresolveLinTerm[(1.0, 1)]),
        -2.0,
        5.0,
    )

    @test PC.get_quad_coeff(con.qe, 1, 2) == 3.0
    @test !PC.is_integer(con)

    PC._scale_constraint_by_two!(con)

    @test PC.get_quad_coeff(con.qe, 1, 2) == 6.0
    @test PC.get_lin_coeff(con.qe, 1) == 2.0
    @test con.lhs == -4.0
    @test con.rhs == 10.0
    @test PC.is_integer(con)
end

@testset "presolve! handles odd bilinear coefficients" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(0.0, 1.0),
        2 => PC.IntVar(0.0, 1.0),
    )
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[(3.0, 1, 2)], PresolveLinTerm[]),
        3.0,
        3.0,
    )
    model = PC.QPModel(vars, [con], presolve_empty_objective(), :min)

    @test !PC.is_integer(con)

    result = QIPresolve.presolve!(
        model;
        residue_strategy = :powers_of_two,
        residue_threshold = 2,
    )

    @test result.model === model
    @test !model.infeasible
    @test all(PC.is_integer, model.cons)
end

@testset "presolve! initially normalizes fixed variables with postsolve data" begin
    vars = Dict{PC.VarId, PC.IntVar}(
        1 => PC.IntVar(2.0, 2.0),
        2 => PC.IntVar(0.0, 5.0),
    )
    con = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[], PresolveLinTerm[(1.0, 1), (1.0, 2)]),
        7.0,
        7.0,
    )
    con.qe.constant = 2.0
    model = PC.QPModel(vars, [con], presolve_empty_objective(), :min)

    result = QIPresolve.presolve!(
        model;
        residue_strategy = :powers_of_two,
        residue_threshold = 2,
    )

    @test !model.infeasible
    @test isempty(model.vars)
    @test isempty(model.cons)
    @test QIPresolve.postsolve(result.postsolver, Dict{PC.VarId, Float64}()) ==
        Dict{PC.VarId, Float64}(1 => 2.0, 2 => 3.0)
end

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

@testset "presolve! passes parity strategy option" begin
    function conjunction_true_model()
        vars = Dict{PC.VarId, PC.IntVar}(
            1 => PC.IntVar(0.0, 1.0),
            2 => PC.IntVar(0.0, 1.0),
        )
        con = PC.Constraint(
            presolve_next_con_id(),
            PC.QuadExpr(PresolveQuadTerm[(2.0, 1, 2)], PresolveLinTerm[]),
            2.0,
            2.0,
        )
        return PC.QPModel(vars, [con], presolve_empty_objective(), :min)
    end

    mod2_model = conjunction_true_model()
    mod2_result = QIPresolve.presolve!(
        mod2_model;
        parity_strategy = :mod2_basic,
        residue_strategy = :powers_of_two,
        residue_threshold = 0,
    )

    @test mod2_result.model === mod2_model
    @test !mod2_model.infeasible
    @test length(mod2_model.vars) == 2

    full_model = conjunction_true_model()
    full_result = QIPresolve.presolve!(
        full_model;
        parity_strategy = :full,
        residue_strategy = :powers_of_two,
        residue_threshold = 0,
    )

    @test full_result.model === full_model
    @test !full_model.infeasible
    @test isempty(full_model.vars)

    @test_throws ArgumentError QIPresolve.presolve!(
        conjunction_true_model();
        parity_strategy = :invalid,
    )
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
        1 => PC.IntVar(0.0, 4.0),
        2 => PC.IntVar(0.0, 4.0),
        3 => PC.IntVar(0.0, 4.0),
        4 => PC.IntVar(0.0, 4.0),
    )
    con1 = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[(2.0, 1, 3)], PresolveLinTerm[]),
        1.0,
        9.0,
    )
    con2 = PC.Constraint(
        presolve_next_con_id(),
        PC.QuadExpr(PresolveQuadTerm[(2.0, 2, 4)], PresolveLinTerm[]),
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

    stats = PC.parity_presolve!(model, propagator, postsolver)
    new_con_ids = setdiff(Set(con.id for con in model.cons), Set([con.id]))

    @test stats.changed
    @test stats.domains_changed
    @test stats.pattern_rewritten_vars == 2
    @test !isempty(new_con_ids)
    @test all(con_id in stats.coefficient_changed_constraint_ids for con_id in new_con_ids)
end

@testset "presolve! handles stale parity events from script random model" begin
    random_qip_kwargs = (
        p_con_eq = 0.0,
        var_threshold_lb = -10,
        var_threshold_ub = 10,
        p_var_is_candidate = 0.02,
        p_var_bilin = 0.4,
        p_var_diag = 0.5,
        p_var_lin = 0.0,
        coeff_lb = -50,
        coeff_ub = 50,
        force_diag_even = false,
        force_lin_even = false,
        force_feasibility = true,
        constraint_slack_range = -10:10,
    )
    jump_model, _ = generate_random_qip_model(100, 200; random_qip_kwargs..., seed = 42)
    model = QIPresolve.build_model(QIPresolve.from_moi(backend(jump_model)))

    result = QIPresolve.presolve!(model)

    @test result isa QIPresolve.PresolveResult
    @test result.model === model
    @test !model.infeasible
end
