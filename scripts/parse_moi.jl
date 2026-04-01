import QIPresolve as QIP

moi = QIP.load_moi_model("small.mof.json")


qp_model_builder = QIP.from_moi(moi)

qp_model = QIP.build_model(qp_model_builder)

QIP.fix_vars!(qp_model)


xor_model = QIP.PresolvingCore.build_parity_model(qp_model)


QIP.PresolvingCore.reformulate_bipartite_cons!(xor_model)


QIP.PresolvingCore.propagate!(xor_model)

for con in xor_model.cons
    if con.is_pure_xor
        println(con)
    end
end

for con in xor_model.cons
    if !con.is_pure_xor
        println(con)
    end
end

QIP.PresolvingCore.gauss_jordan!(xor_model; skip_conj = true)


#look at waht happens to c9?? var removal does not work like expected
