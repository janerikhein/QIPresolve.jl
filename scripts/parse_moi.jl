using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

import QIPresolve as QIP

PC = QIP.PresolvingCore

moi = QIP.load_moi_model("large2.mof.json")
qp_model_builder = QIP.from_moi(moi)
qp_model = QIP.build_model(qp_model_builder)
PC.fix_vars!(qp_model)

log_domain_sum(model::PC.QPModel) = sum(log(var.ub - var.lb + 1) for var in values(model.vars))

println("nvars=$(length(qp_model.vars)) ncons=$(length(qp_model.cons))")
println("log_domain_sum_before = $(log_domain_sum(qp_model))")

propagator = PC.PropagationManager(Int[])
phase_idx = 0

println(qp_model)

while !qp_model.infeasible
    global phase_idx += 1
    stats = PC.parity_presolve_phase!(qp_model, propagator)
    PC.scale_constraints_gcd!(qp_model)
    println(
        "phase $phase_idx: changed=$(stats.changed) " *
        "fixed_parities=$(stats.fixed_parities) " *
        "pattern_rewritten_vars=$(stats.pattern_rewritten_vars) " *
        "infeasible=$(qp_model.infeasible) " *
        "nvars=$(length(qp_model.vars)) ncons=$(length(qp_model.cons))"
    )

    stats.changed || break
end

println(qp_model)
println("log_domain_sum_after = $(log_domain_sum(qp_model))")
println("final model infeasible = $(qp_model.infeasible)")
