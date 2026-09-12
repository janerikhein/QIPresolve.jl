using JuMP, SCIP

length(ARGS) == 1 || error("Usage: julia --project=. scripts/solve_lp_scip.jl FILE.lp")
isfile(ARGS[1]) || error("File not found: $(ARGS[1])")

model = read_from_file(ARGS[1])
set_optimizer(model, SCIP.Optimizer)
optimize!(model)
println("Termination status: ", termination_status(model))
