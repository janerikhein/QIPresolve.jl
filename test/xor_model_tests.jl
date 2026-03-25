using Test
using Random
import QIPresolve.PresolvingCore as PC

const XorModelQuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const XorModelLinTerm = Tuple{Float64, PC.VarId}

function empty_objective()
    return PC.QuadExpr(XorModelQuadTerm[], XorModelLinTerm[]; constant = 0.0)
end

function count_conj_terms(con::PC.XorConstraint)
    con.conj === nothing && return 0
    n = size(con.conj, 1)
    total = 0
    for i in 1:n
        for j in (i + 1):n
            con.conj[i, j] && (total += 1)
        end
    end
    return total
end

function parity_constraint_holds(con::PC.XorConstraint, pos_to_var::Vector{PC.VarId}, x::Vector{Int})
    lhs = false
    n = length(con.par)

    for pos in 1:n
        con.par[pos] || continue
        lhs = xor(lhs, isodd(x[pos_to_var[pos]]))
    end

    if con.conj !== nothing
        for i in 1:n
            xi_odd = isodd(x[pos_to_var[i]])
            for j in (i + 1):n
                con.conj[i, j] || continue
                lhs = xor(lhs, xi_odd && isodd(x[pos_to_var[j]]))
            end
        end
    end

    return lhs == con.rhs
end

function all_parity_constraints_hold(model::PC.ParityModel, x::Vector{Int})
    for con in model.cons
        parity_constraint_holds(con, model.pos_to_var_id, x) || return false
    end
    return true
end

function refresh_constraint_cache!(model::PC.ParityModel)
    for con in model.cons
        PC.update!(con)
    end
    return model
end

function random_equality_model_with_solution(rng::AbstractRNG, nvars::Int, ncons::Int)
    vars = Dict{PC.VarId, PC.IntVar}()
    x = zeros(Int, nvars)

    for vid in 1:nvars
        if rand(rng) < 0.5
            vars[vid] = PC.IntVar(0.0, 1.0)
            x[vid] = rand(rng, 0:1)
        else
            lb = rand(rng, -3:0)
            ub = rand(rng, 1:4)
            vars[vid] = PC.IntVar(lb, ub)
            x[vid] = rand(rng, lb:ub)
        end
    end

    cons = Vector{PC.Constraint}(undef, ncons)
    for cid in 1:ncons
        quad_terms = XorModelQuadTerm[]
        lin_terms = XorModelLinTerm[]

        for i in 1:nvars
            diag = rand(rng, -6:6)
            diag == 0 || push!(quad_terms, (diag, i, i))

            lin = rand(rng, -6:6)
            lin == 0 || push!(lin_terms, (lin, i))

            for j in (i + 1):nvars
                qij = rand(rng, -8:2:8)
                qij == 0 || push!(quad_terms, (qij, i, j))
            end
        end

        if isempty(quad_terms) && isempty(lin_terms)
            vid = rand(rng, 1:nvars)
            push!(lin_terms, (1.0, vid))
        end

        qe = PC.QuadExpr(quad_terms, lin_terms; constant = 0.0)
        rhs = PC.eval_full(qe, x)
        cons[cid] = PC.Constraint(cid, qe, rhs, rhs)
    end

    model = PC.QPModel(vars, cons, empty_objective(), :min)
    return model, x
end

@testset "get_parity_model basic extraction" begin
    vars = Dict{PC.VarId, PC.IntVar}(1 => PC.IntVar(0.0, 1.0), 2 => PC.IntVar(0.0, 1.0))
    qe = PC.QuadExpr(XorModelQuadTerm[], XorModelLinTerm[(1.0, 1), (1.0, 2)])
    con = PC.Constraint(1, qe, 1.0, 1.0)
    model = PC.QPModel(vars, [con], empty_objective(), :min)

    parity_model = PC.get_parity_model(model)
    @test length(parity_model.cons) == 2
    @test Set(parity_model.pos_to_var_id) == Set([1, 2])

    pure = [c for c in parity_model.cons if c.meta.is_pure_xor]
    mixed = [c for c in parity_model.cons if !c.meta.is_pure_xor]
    @test length(pure) == 1
    @test length(mixed) == 1

    @test count(pure[1].par) == 2
    @test pure[1].rhs
    @test count(mixed[1].par) == 0
    @test !mixed[1].rhs
    @test count_conj_terms(mixed[1]) == 1

    x = [1, 0]
    @test all_parity_constraints_hold(parity_model, x)
end

@testset "propagate! splits full-bipartite conjunction rows" begin
    pos_to_var = [1, 2, 3, 4]
    var_to_pos = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var))
    conj = falses(4, 4)
    for (i, j) in ((1, 3), (1, 4), (2, 3), (2, 4))
        conj[i, j] = true
        conj[j, i] = true
    end

    con = PC.XorConstraint(falses(4), conj, true)
    model = PC.ParityModel(var_to_pos, pos_to_var, [con])

    x = [1, 0, 1, 0]
    @test all_parity_constraints_hold(model, x)

    status = PC.propagate!(model)
    @test status == PC.PARITY_PROPAGATE_UPDATED
    @test length(model.cons) == 2
    @test all(c -> c.meta.is_pure_xor, model.cons)
    @test all_parity_constraints_hold(model, x)
end

@testset "gauss_jordan! preserves validity for known solutions" begin
    pos_to_var = [1, 2, 3]
    var_to_pos = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var))
    cons = PC.XorConstraint[
        PC.XorConstraint(BitVector([1, 1, 1]), true),
        PC.XorConstraint(BitVector([1, 1, 0]), false),
        PC.XorConstraint(BitVector([0, 1, 1]), true),
    ]
    model = PC.ParityModel(var_to_pos, pos_to_var, cons)
    x = [0, 0, 1]

    @test all_parity_constraints_hold(model, x)
    PC.gauss_jordan!(model)
    @test all_parity_constraints_hold(model, x)

    pos_to_var2 = [1, 2]
    var_to_pos2 = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var2))
    conj = falses(2, 2)
    conj[1, 2] = true
    conj[2, 1] = true
    model2 = PC.ParityModel(var_to_pos2, pos_to_var2, PC.XorConstraint[
        PC.XorConstraint(BitVector([1, 0]), true),
        PC.XorConstraint(falses(2), conj, false),
    ])
    x2 = [1, 0]

    @test all_parity_constraints_hold(model2, x2)
    PC.gauss_jordan!(model2; skip_conj = true)
    @test all_parity_constraints_hold(model2, x2)
end

@testset "random equality QPModels keep parity validity invariant" begin
    rng = MersenneTwister(20260301)

    for _ in 1:60
        nvars = rand(rng, 2:7)
        ncons = rand(rng, 1:10)
        qp_model, x = random_equality_model_with_solution(rng, nvars, ncons)

        parity_model = PC.get_parity_model(qp_model)
        @test all_parity_constraints_hold(parity_model, x)

        model_after_propagation = deepcopy(parity_model)
        PC.propagate!(model_after_propagation)
        @test all_parity_constraints_hold(model_after_propagation, x)

        model_after_elimination = deepcopy(model_after_propagation)
        refresh_constraint_cache!(model_after_elimination)
        PC.gauss_jordan!(model_after_elimination)
        @test all_parity_constraints_hold(model_after_elimination, x)

        model_skip_conj = deepcopy(parity_model)
        refresh_constraint_cache!(model_skip_conj)
        PC.gauss_jordan!(model_skip_conj; skip_conj = true)
        @test all_parity_constraints_hold(model_skip_conj, x)
    end
end
