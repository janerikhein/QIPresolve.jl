# Included inside MainBenchmarkPresolveComparison. All native pointers stay
# inside run_scip; captured constraints and the optimizer are released in finally.

const NativeVar = Ptr{SCIP.SCIP_VAR}
const NativeCons = Ptr{SCIP.SCIP_CONS}
const Polynomial = Dict{Tuple{UInt, UInt}, Float64}
const CONSTANT = (UInt(0), UInt(0))

function add_polynomial!(target::Polynomial, source::Polynomial, scale::Float64 = 1.0)
    for (key, value) in source
        target[key] = get(target, key, 0.0) + scale * value
    end
    filter!(p -> !iszero(last(p)), target)
    return target
end

function multiply_polynomials(left::Polynomial, right::Polynomial)
    result = Polynomial()
    for (a, av) in left, (b, bv) in right
        variables = sort!(filter(!iszero, [a..., b...]))
        length(variables) <= 2 || error("SCIP expression is no longer quadratic")
        key = isempty(variables) ? CONSTANT : length(variables) == 1 ?
            (variables[1], UInt(0)) : (variables[1], variables[2])
        result[key] = get(result, key, 0.0) + av * bv
    end
    return result
end

function expression_polynomial(expr)
    handler = unsafe_string(SCIP.SCIPexprhdlrGetName(SCIP.SCIPexprGetHdlr(expr)))
    if handler == "var"
        return Polynomial((UInt(SCIP.SCIPgetVarExprVar(expr)), UInt(0)) => 1.0)
    elseif handler == "val"
        return Polynomial(CONSTANT => SCIP.SCIPgetValueExprValue(expr))
    end
    children = unsafe_wrap(Array, SCIP.SCIPexprGetChildren(expr), SCIP.SCIPexprGetNChildren(expr))
    if handler == "sum"
        result = Polynomial(CONSTANT => SCIP.SCIPgetConstantExprSum(expr))
        coefs = unsafe_wrap(Array, SCIP.SCIPgetCoefsExprSum(expr), length(children))
        for (child, coef) in zip(children, coefs)
            add_polynomial!(result, expression_polynomial(child), coef)
        end
        return result
    elseif handler == "prod"
        result = Polynomial(CONSTANT => SCIP.SCIPgetCoefExprProduct(expr))
        for child in children
            result = multiply_polynomials(result, expression_polynomial(child))
        end
        return result
    elseif handler == "pow"
        exponent = SCIP.SCIPgetExponentExprPow(expr)
        exponent in (0.0, 1.0, 2.0) || error("Unsupported SCIP power: $exponent")
        result = Polynomial(CONSTANT => 1.0)
        for _ in 1:Int(exponent)
            result = multiply_polynomials(result, expression_polynomial(only(children)))
        end
        return result
    end
    error("Unsupported SCIP expression handler: $handler")
end

function linear_polynomial(variables, coefficients)
    result = Polynomial()
    for (var, coef) in zip(variables, coefficients)
        key = (UInt(var), UInt(0))
        result[key] = get(result, key, 0.0) + Float64(coef)
    end
    return result
end

function native_constraint(scip, con)
    handler = unsafe_string(SCIP.SCIPconshdlrGetName(SCIP.SCIPconsGetHdlr(con)))
    if handler == "nonlinear"
        return (poly = expression_polynomial(SCIP.SCIPgetExprNonlinear(con)),
            lhs = SCIP.SCIPgetLhsNonlinear(con), rhs = SCIP.SCIPgetRhsNonlinear(con))
    elseif handler == "linear"
        n = SCIP.SCIPgetNVarsLinear(scip, con)
        vars = unsafe_wrap(Array, SCIP.SCIPgetVarsLinear(scip, con), n)
        vals = unsafe_wrap(Array, SCIP.SCIPgetValsLinear(scip, con), n)
        return (poly = linear_polynomial(vars, vals), lhs = SCIP.SCIPgetLhsLinear(scip, con),
            rhs = SCIP.SCIPgetRhsLinear(scip, con))
    elseif handler == "varbound"
        vars = [SCIP.SCIPgetVarVarbound(scip, con), SCIP.SCIPgetVbdvarVarbound(scip, con)]
        return (poly = linear_polynomial(vars, [1.0, SCIP.SCIPgetVbdcoefVarbound(scip, con)]),
            lhs = SCIP.SCIPgetLhsVarbound(scip, con), rhs = SCIP.SCIPgetRhsVarbound(scip, con))
    elseif handler == "knapsack"
        n = SCIP.SCIPgetNVarsKnapsack(scip, con)
        vars = unsafe_wrap(Array, SCIP.SCIPgetVarsKnapsack(scip, con), n)
        weights = unsafe_wrap(Array, SCIP.SCIPgetWeightsKnapsack(scip, con), n)
        return (poly = linear_polynomial(vars, weights), lhs = -Inf,
            rhs = Float64(SCIP.SCIPgetCapacityKnapsack(scip, con)))
    elseif handler == "setppc"
        n = SCIP.SCIPgetNVarsSetppc(scip, con)
        vars = unsafe_wrap(Array, SCIP.SCIPgetVarsSetppc(scip, con), n)
        kind = SCIP.SCIPgetTypeSetppc(scip, con)
        lhs = kind == SCIP.SCIP_SETPPCTYPE_PACKING ? -Inf : 1.0
        rhs = kind == SCIP.SCIP_SETPPCTYPE_COVERING ? Inf : 1.0
        return (poly = linear_polynomial(vars, ones(n)), lhs = lhs, rhs = rhs)
    elseif handler == "logicor"
        n = SCIP.SCIPgetNVarsLogicor(scip, con)
        vars = unsafe_wrap(Array, SCIP.SCIPgetVarsLogicor(scip, con), n)
        return (poly = linear_polynomial(vars, ones(n)), lhs = 1.0, rhs = Inf)
    elseif handler == "and"
        n = SCIP.SCIPgetNVarsAnd(scip, con)
        vars = unsafe_wrap(Array, SCIP.SCIPgetVarsAnd(scip, con), n)
        poly = Polynomial(CONSTANT => 1.0)
        for var in vars
            poly = multiply_polynomials(poly, linear_polynomial([var], [1.0]))
        end
        add_polynomial!(poly, linear_polynomial([SCIP.SCIPgetResultantAnd(scip, con)], [-1.0]))
        return (poly = poly, lhs = 0.0, rhs = 0.0)
    end
    error("Unsupported SCIP constraint handler: $handler")
end

"Resolve fixed, negated, and (multi-)aggregated variables into active variables."
function active_variable_polynomial(scip, variable::NativeVar)
    size = 8
    while true
        vars = Vector{NativeVar}(undef, size)
        coefficients = zeros(size)
        vars[1] = variable
        coefficients[1] = 1.0
        n = Ref{Cint}(1)
        constant = Ref{Cdouble}(0.0)
        required = Ref{Cint}(0)
        SCIP.@SCIP_CALL SCIP.SCIPgetProbvarLinearSum(scip, vars, coefficients, n,
            size, constant, required)
        if required[] > size
            size = Int(required[])
            continue
        end
        poly = linear_polynomial(view(vars, 1:n[]), view(coefficients, 1:n[]))
        poly[CONSTANT] = constant[]
        return poly
    end
end

function active_polynomial(scip, poly::Polynomial, cache)
    result = Polynomial()
    for ((a, b), coefficient) in poly
        term = Polynomial(CONSTANT => coefficient)
        for address in (a, b)
            address == 0 && continue
            substitute = get!(cache, address) do
                active_variable_polynomial(scip, NativeVar(address))
            end
            term = multiply_polynomials(term, substitute)
        end
        add_polynomial!(result, term)
    end
    # SCIP may replace binary squares by their linear equivalent.
    for ((a, b), value) in collect(result)
        if a != 0 && a == b && SCIP.SCIPvarGetType(NativeVar(a)) == SCIP.SCIP_VARTYPE_BINARY
            delete!(result, (a, b))
            result[(a, UInt(0))] = get(result, (a, UInt(0)), 0.0) + value
        end
    end
    return result
end

"Return s when after = s * before + constant; shifts do not change widths."
function proportional_scale(before::Polynomial, after::Polynomial)
    keys_union = union(keys(before), keys(after))
    delete!(keys_union, CONSTANT)
    isempty(keys_union) && return nothing
    key = argmax(k -> abs(get(before, k, 0.0)), collect(keys_union))
    coefficient = get(before, key, 0.0)
    coefficient == 0.0 && return nothing
    scale = get(after, key, 0.0) / coefficient
    isfinite(scale) && scale != 0.0 || return nothing
    for k in keys_union
        isapprox(get(after, k, 0.0), scale * get(before, k, 0.0); rtol = 1e-8, atol = 1e-8) || return nothing
    end
    return scale
end

native_name(con) = unsafe_string(SCIP.SCIPconsGetName(con))

function build_native_model(model)
    optimizer = SCIP.Optimizer()
    try
        variable_map = QIP.ModelIO._register_moi_variables!(optimizer, model.vars)
        con_map = IdDict{PC.Constraint, NativeCons}()
        for (index, con) in enumerate(model.cons)
            func = QIP.ModelIO._moi_scalar_function(con.qe, variable_map)
            set = QIP.ModelIO._moi_bound_set(con.lhs, con.rhs)
            set === nothing && continue
            ci = MOI.add_constraint(optimizer, func, set)
            MOI.set(optimizer, MOI.ConstraintName(), ci, "qipc_$index")
            con_map[con] = SCIP.cons(optimizer, ci)
        end
        obj = QIP.ModelIO._moi_scalar_function(model.obj_expr, variable_map)
        MOI.set(optimizer, MOI.ObjectiveSense(), QIP.ModelIO.QPObjSenseMapping[model.obj_sense])
        MOI.set(optimizer, MOI.ObjectiveFunction{typeof(obj)}(), obj)
        return optimizer, con_map
    catch
        SCIP.free_scip(optimizer.inner)
        rethrow()
    end
end

function run_scip(model, baselines, config; diagnostics::IO = stderr, label = "scip")
    contributions = core_contributions(model, baselines)
    if model.infeasible
        return (status = "infeasible", log_domain = 0.0, contributions = contributions, nodes = 0)
    end
    optimizer, con_map = build_native_model(model)
    scip = optimizer.inner
    captured = NativeCons[]
    try
        config.scip_config === nothing || SCIP.@SCIP_CALL SCIP.SCIPreadParams(scip, config.scip_config)
        MOI.set(optimizer, MOI.Silent(), true)
        SCIP.@SCIP_CALL SCIP.SCIPtransformProb(scip)
        tracked = []
        for (index, baseline) in enumerate(baselines)
            haskey(con_map, baseline.con) || continue
            transformed = Ref{NativeCons}()
            SCIP.@SCIP_CALL SCIP.SCIPgetTransformedCons(scip, con_map[baseline.con], transformed)
            SCIP.@SCIP_CALL SCIP.SCIPcaptureCons(scip, transformed[])
            push!(captured, transformed[])
            push!(tracked, (index = index, con = transformed[], name = native_name(transformed[]),
                before = native_constraint(scip, transformed[]),
                core_scale = baseline.scale / baseline.con._bound_scale))
        end
        source_names = Set(native_name(c) for c in values(con_map))
        SCIP.@SCIP_CALL SCIP.SCIPpresolve(scip)
        nodes = Int(SCIP.SCIPgetNNodes(scip))
        nodes == 0 || error("SCIPpresolve unexpectedly processed $nodes search nodes")
        native_status = SCIP.SCIPgetStatus(scip)
        active = [c for c in unsafe_wrap(Array, SCIP.SCIPgetConss(scip), SCIP.SCIPgetNConss(scip))
            if SCIP.SCIPconsIsDeleted(c) == 0]
        active_set = Set(active)
        names = Dict{String, Vector{NativeCons}}()
        for con in active
            push!(get!(names, native_name(con), NativeCons[]), con)
        end
        # Reserve surviving source constraints before matching renamed upgrades;
        # an aggregated-away duplicate must not claim another input row,
        # including helper constraints introduced by QIPresolve.
        reserved = Set(c for c in active if native_name(c) in source_names || any(t -> t.con == c, tracked))
        claimed = Set{NativeCons}()
        cache = Dict{UInt, Polynomial}()
        after_data = Dict{NativeCons, Any}()
        for con in active
            after_data[con] = try
                data = native_constraint(scip, con)
                merge(data, (active_poly = active_polynomial(scip, data.poly, cache),))
            catch err
                println(diagnostics, "$label: cannot inspect $(native_name(con)): $(sprint(showerror, err))")
                nothing
            end
        end
        for track in tracked
            before_poly = active_polynomial(scip, track.before.poly, cache)
            candidates = track.con in active_set ? [track.con] : get(names, track.name, NativeCons[])
            if isempty(candidates)
                # Upgrades normally preserve names. A changed name is accepted
                # only with a verified proportional expression and unique owner.
                candidates = [c for c in active if !(c in reserved) && !(c in claimed) &&
                    after_data[c] !== nothing &&
                    proportional_scale(before_poly, after_data[c].active_poly) !== nothing]
            end
            if isempty(candidates)
                if SCIP.SCIPconsIsDeleted(track.con) != 0
                    contributions[track.index] = 1.0
                else
                    println(diagnostics, "$label: untraceable $(track.name); retaining last verified interval")
                end
                continue
            end
            baseline = baselines[track.index]
            if !isfinite(baseline.lhs) || !isfinite(baseline.rhs) || baseline.lhs == baseline.rhs
                # Surviving equalities/one-sided rows always contribute zero;
                # their coefficient rewrites cannot affect this metric.
                union!(claimed, candidates)
                contributions[track.index] = 0.0
                continue
            end
            # If a range is split during an upgrade, intersect its successor rows
            # in the original expression's units, including affine offsets.
            lhs, rhs = track.before.lhs, track.before.rhs
            verified = false
            for candidate in candidates
                push!(claimed, candidate)
                after = after_data[candidate]
                after === nothing && continue
                scale = proportional_scale(before_poly, after.active_poly)
                scale === nothing && continue
                offset = get(after.active_poly, CONSTANT, 0.0) - scale * get(before_poly, CONSTANT, 0.0)
                lo, hi = (after.lhs - offset) / scale, (after.rhs - offset) / scale
                scale < 0 && ((lo, hi) = (hi, lo))
                lhs, rhs = max(lhs, lo), min(rhs, hi)
                verified = true
            end
            if verified
                contributions[track.index] = bound_contribution(baselines[track.index],
                    (rhs - lhs) * track.core_scale)
            else
                println(diagnostics, "$label: nonproportional rewrite of $(track.name); retaining last verified interval")
            end
        end
        infeasible = native_status == SCIP.SCIP_STATUS_INFEASIBLE
        vars = unsafe_wrap(Array, SCIP.SCIPgetVars(scip), SCIP.SCIPgetNVars(scip))
        log_domain = infeasible ? 0.0 : sum(vars; init = 0.0) do var
            lb, ub = SCIP.SCIPvarGetLbGlobal(var), SCIP.SCIPvarGetUbGlobal(var)
            SCIP.SCIPisInfinity(scip, abs(lb)) == 0 && SCIP.SCIPisInfinity(scip, abs(ub)) == 0 ||
                error("Unbounded active SCIP variable: domain metric is undefined")
            lb <= ub || error("Inconsistent SCIP bounds without an infeasibility status")
            log(ub - lb + 1.0)
        end
        status = infeasible ? "infeasible" :
            (native_status == SCIP.SCIP_STATUS_OPTIMAL || isempty(active)) ? "feasible" : "reduced"
        println(diagnostics, "$label: SCIP status=$native_status, search_nodes=$nodes")
        return (status = status, log_domain = log_domain, contributions = contributions, nodes = nodes)
    finally
        for con in captured
            SCIP.@SCIP_CALL SCIP.SCIPreleaseCons(scip, Ref(con))
        end
        SCIP.free_scip(scip)
    end
end
