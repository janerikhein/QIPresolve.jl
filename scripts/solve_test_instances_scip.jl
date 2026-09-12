using JuMP, SCIP, Printf
import MathOptInterface as MOI

const TIME_LIMIT = 60 * 60.0
const INSTANCE_DIR = joinpath(@__DIR__, "..", "instances", "test")

files = filter(endswith(".lp"), readdir(INSTANCE_DIR; join = true))
sort!(files; by = file -> parse(Int, match(r"_n(\d+)\.lp$", file).captures[1]))

runs = Tuple{String, Float64, MOI.TerminationStatusCode}[]
for file in files
    println("\n", "="^80, "\nSolving ", basename(file), "\n", "="^80)
    model = read_from_file(file)
    set_optimizer(model, SCIP.Optimizer)
    set_time_limit_sec(model, TIME_LIMIT)
    optimize!(model)

    status = termination_status(model)
    push!(runs, (basename(file), solve_time(model), status))
    status == MOI.TIME_LIMIT && break
end

println("\nRuntime summary")
@printf("%-24s %12s  %s\n", "instance", "seconds", "status")
for (instance, runtime, status) in runs
    @printf("%-24s %12.3f  %s\n", instance, runtime, string(status))
end
