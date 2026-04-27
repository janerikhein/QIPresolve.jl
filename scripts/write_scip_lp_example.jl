using QIPresolve
import QIPresolve.PresolvingCore as PC
import SCIP

const QuadTerm = Tuple{Float64, PC.VarId, PC.VarId}
const LinTerm = Tuple{Float64, PC.VarId}

vars = Dict{PC.VarId, PC.IntVar}(var_id => PC.IntVar(0.0, 1.0) for var_id in 1:24)

constraints = PC.Constraint[]
constraint_id = 1

for var_id in 1:23
    push!(
        constraints,
        PC.Constraint(
            constraint_id,
            PC.QuadExpr(QuadTerm[], LinTerm[(1.0, var_id), (1.0, var_id + 1)]),
            -Inf,
            1.0,
        ),
    )
    global constraint_id += 1
end

for var_id in 1:4:17
    push!(
        constraints,
        PC.Constraint(
            constraint_id,
            PC.QuadExpr(
                QuadTerm[(1.0, var_id, var_id + 2)],
                LinTerm[(1.0, var_id + 4)],
            ),
            -Inf,
            1.0,
        ),
    )
    global constraint_id += 1
end

objective = PC.QuadExpr(
    QuadTerm[],
    LinTerm[(Float64((var_id % 5) + 1), var_id) for var_id in 1:24];
    constant = 0.0,
)
model = PC.QPModel(vars, constraints, objective, :max)

moi_model = QIPresolve.build_moi_model(model, :scip)

original_path = joinpath("results", "scip_example_original.lp")
presolved_path = joinpath("results", "scip_example_presolved.lp")
mkpath(dirname(original_path))

QIPresolve.ModelIO.save_moi(moi_model, original_path)
println("Wrote $(original_path) ($(filesize(original_path)) bytes)")

backend = getfield(moi_model, :model)
scip = backend.inner

SCIP.@SCIP_CALL SCIP.SCIPpresolve(scip)
SCIP.@SCIP_CALL SCIP.SCIPwriteTransProblem(scip, presolved_path, "lp", true)
println("Wrote $(presolved_path) ($(filesize(presolved_path)) bytes)")
