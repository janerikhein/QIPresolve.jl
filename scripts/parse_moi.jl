using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

import QIPresolve as QIP

PC = QIP.PresolvingCore

moi = QIP.load_moi_model("large.mof.json")


qp_model_builder = QIP.from_moi(moi)

qp_model = QIP.build_model(qp_model_builder)

QIP.fix_vars!(qp_model)


xor_model = PC.build_parity_model(qp_model)


prop = PC.PropagationManager(collect(keys(xor_model.var_id_to_pos)))

PC.reformulate_bipartite_cons!(xor_model)


QIP.PresolvingCore.propagate!(xor_model, prop)



#PC.gauss_jordan_xor!(xor_model)

#PC.substitute_pivots_in_conjunctive_terms!(xor_model)
#PC.cleanup!(xor_model)
#PC.substitute_parity_pivots!(xor_model)


#PC.gauss_jordan_xor_and!(xor_model)
#changed = falses(length(xor_model.cons))
#begin
#    substitution = PC.pop_substitution!(prop)
#    #substitution === nothing && break
#    vid, substid, neg = substitution
#end
#PC._substitute_var_rows!(changed, xor_model, vid, substid, neg)



