# A forbidden assignment (a, b) gives the clause ¬a ∨ ¬b. Using the same
# literal twice gives a unit clause. Register only variables used by a clause.
function _forbid_literals!(manager::PropagationManager, a::VarLit, b::VarLit = a)
    ensure_literals!(manager, (a.vid, b.vid))
    add_implication!(manager, a, negated(b))
end

function _seed_binary_pair!(manager, con, x, y)
    q, a = quad(con.qe), lin(con.qe)
    all(isfinite, q) && all(isfinite, a) || return
    # Use exact values of the stored Float64 coefficients, including both
    # off-diagonal entries as in eval_full. Rounded sums can exclude a corner.
    exact = Rational{BigInt}
    c = exact(con.qe.constant)
    ax = exact(a[1]) + exact(q[1, 1])
    ay = exact(a[2]) + exact(q[2, 2])
    axy = exact(q[1, 2]) + exact(q[2, 1])
    lo, hi = exact(con.lhs), exact(con.rhs)
    allowed = map(v -> lo <= v <= hi, (c, c + ay, c + ax, c + ax + ay + axy))
    all(allowed) && return
    if !any(allowed)
        manager.infeasible = true
        return
    end
    f00, f01, f10, f11 = allowed
    x0, x1 = f00 || f01, f10 || f11
    y0, y1 = f00 || f10, f01 || f11
    if !(x0 && x1 && y0 && y1)
        x0 || _forbid_literals!(manager, _lit_for_value(x, false))
        x1 || _forbid_literals!(manager, _lit_for_value(x, true))
        y0 || _forbid_literals!(manager, _lit_for_value(y, false))
        y1 || _forbid_literals!(manager, _lit_for_value(y, true))
        return
    end
    for (i, feasible) in enumerate(allowed)
        feasible && continue
        _forbid_literals!(manager, _lit_for_value(x, i > 2), _lit_for_value(y, iseven(i)))
    end
end

# Terms are sorted by decreasing cost. Once a pair fits the slack, all
# remaining partners fit too. A unit exclusion subsumes every pair containing it.
function _seed_binary_bound!(manager, terms, slack, lower::Bool)
    if slack < 0
        manager.infeasible = true
        return
    end
    adverse(lit) = lower ? negated(lit) : lit
    for (i, (cost, lit)) in enumerate(terms)
        2 * cost <= slack && break
        if cost > slack
            _forbid_literals!(manager, adverse(lit))
            continue
        end
        for j in (i + 1):length(terms)
            other_cost, other = terms[j]
            cost + other_cost <= slack && break
            _forbid_literals!(manager, adverse(lit), adverse(other))
        end
    end
end

"""Seed unit and binary clauses implied by the bounds of normalized all-binary rows."""
function _seed_binary_implications!(manager::PropagationManager, model::QPModel)
    for con in model.cons
        ids = vars(con)
        length(ids) >= 2 && all(id -> is_binary(model.vars[id]), ids) || continue
        isfinite(con.qe.constant) && !isnan(con.lhs) && !isnan(con.rhs) || continue
        if con.lhs > con.rhs
            manager.infeasible = true
            break
        end
        if length(ids) == 2
            _seed_binary_pair!(manager, con, ids...)
        elseif all(iszero, quad(con.qe)) && all(isfinite, lin(con.qe))
            minimum_activity = maximum_activity = Rational{BigInt}(con.qe.constant)
            terms = Tuple{Rational{BigInt}, VarLit}[]
            for (vid, coefficient) in zip(ids, lin(con.qe))
                a = Rational{BigInt}(coefficient)
                minimum_activity += min(a, 0)
                maximum_activity += max(a, 0)
                push!(terms, (abs(a), _lit_for_value(vid, a > 0)))
            end
            sort!(terms; by = first, rev = true)
            _seed_binary_bound!(manager, terms, Rational{BigInt}(con.rhs) - minimum_activity, false)
            _seed_binary_bound!(manager, terms, maximum_activity - Rational{BigInt}(con.lhs), true)
        end
        manager.infeasible && break
    end
    return manager
end
