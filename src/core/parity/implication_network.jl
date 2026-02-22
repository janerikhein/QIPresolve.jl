using Graphs

"""
    ImpGraph(n)

Implication graph for `n` boolean parities p1..pn.
Nodes are literals:
  - 1..n        =>  p1..pn
  - n+1..2n     =>  ¬p1..¬pn

Invariant enforced by `add_implication!`:
  if (X → Y) then (¬Y → ¬X).
"""
mutable struct ImpGraph
    n::Int
    g::DiGraph                 # literal implication graph, size 2n
    dirty::Bool                # SCC/condensation needs recompute?

    # SCC data (valid if !dirty and computed at least once)
    comp_id::Vector{Int}       # comp_id[v] = SCC index of literal node v
    comps::Vector{Vector{Int}} # list of SCCs (nodes)
    cg::DiGraph                # condensation graph (SCC DAG)

    # Reachability cache on cg (optional; computed in recompute!)
    topo::Vector{Int}          # topological order of SCC DAG
    reach::Vector{BitVector}   # reach[c] bitset: SCCs reachable from c (including itself)
end

# ---------------------------
# Literal encoding utilities
# ---------------------------

@inline num_literals(net::ImpGraph) = 2 * net.n

"""
    lit(i; neg=false)

Convenience: literal node id for variable i in 1..n.
If neg=true, returns ¬pi.
"""
@inline function lit(net::ImpGraph, i::Int; neg::Bool=false)
    @boundscheck 1 <= i <= net.n || throw(BoundsError("variable index i=$i out of 1..$(net.n)"))
    return neg ? i + net.n : i
end

"""
    neg(v)

Return node id of negated literal.
"""
@inline function neg(net::ImpGraph, v::Int)
    @boundscheck 1 <= v <= 2net.n || throw(BoundsError("literal v=$v out of 1..$(2net.n)"))
    return v <= net.n ? v + net.n : v - net.n
end

"""
    var_of(v) -> (i, is_neg)

Map literal node id back to variable index and sign.
"""
@inline function var_of(net::ImpGraph, v::Int)
    @boundscheck 1 <= v <= 2net.n || throw(BoundsError("literal v=$v out of 1..$(2net.n)"))
    if v <= net.n
        return (v, false)
    else
        return (v - net.n, true)
    end
end

# ---------------------------
# Construction
# ---------------------------

function ImpGraph(n::Int)
    n >= 1 || throw(ArgumentError("n must be ≥ 1"))
    g = DiGraph(2n)
    return ImpGraph(
        n,
        g,
        true,
        zeros(Int, 2n),
        Vector{Vector{Int}}(),
        DiGraph(0),
        Int[],
        BitVector[],
    )
end

# ---------------------------
# Adding implications
# ---------------------------

"""
    add_implication!(net, x, y)

Add implication x → y, and also add the contrapositive ¬y → ¬x
to enforce the invariant in your slide.

Arguments x,y are literal node ids in 1..2n.
"""
function add_implication!(net::ImpGraph, x::Int, y::Int)
    @boundscheck 1 <= x <= 2net.n || throw(BoundsError("x=$x out of 1..$(2net.n)"))
    @boundscheck 1 <= y <= 2net.n || throw(BoundsError("y=$y out of 1..$(2net.n)"))
    add_edge!(net.g, x, y)
    add_edge!(net.g, neg(net, y), neg(net, x))
    net.dirty = true
    return net
end

"""
    add_equivalence!(net, x, y)

Add x ↔ y (two implications) with invariant.
"""
function add_equivalence!(net::ImpGraph, x::Int, y::Int)
    add_implication!(net, x, y)
    add_implication!(net, y, x)
    return net
end

"""
    add_atmost_one_true!(net, a, b)

Add (a → ¬b) and (b → ¬a), i.e., not (a ∧ b).
Useful for rules like "pa ⇒ ¬pb".
"""
function add_atmost_one_true!(net::ImpGraph, a::Int, b::Int)
    add_implication!(net, a, neg(net, b))
    add_implication!(net, b, neg(net, a))
    return net
end

# ---------------------------
# Recompute SCC + DAG + reachability
# ---------------------------

"""
    recompute!(net; compute_reach=true)

Recompute SCCs, condensation DAG, and optionally a transitive reachability cache
on the SCC DAG using bitsets (fast for repeated path queries).

Call this after a batch of `add_*` operations.
"""
function recompute!(net::ImpGraph; compute_reach::Bool=true)
    if !net.dirty && nv(net.cg) > 0
        return net
    end

    # SCCs of literal graph
    comps = strongly_connected_components(net.g)
    comp_id = zeros(Int, nv(net.g))
    for (cid, comp) in enumerate(comps)
        for v in comp
            comp_id[v] = cid
        end
    end

    # Condensation graph (SCC DAG)
    cg = condensation_graph(net.g)

    # Topological order of SCC DAG
    topo = topological_sort(cg)

    reach = BitVector[]
    if compute_reach
        # Bitset reachability DP in reverse topological order
        k = nv(cg)
        reach = [falses(k) for _ in 1:k]
        for c in reverse(topo)
            r = reach[c]
            r[c] = true
            for nb in outneighbors(cg, c)
                r .|= reach[nb]
            end
        end
    end

    net.comps = comps
    net.comp_id = comp_id
    net.cg = cg
    net.topo = topo
    net.reach = reach
    net.dirty = false
    return net
end

# ---------------------------
# Queries: conflicts, paths, fixings
# ---------------------------

"""
    has_conflict(net) -> Bool

True if ∃i such that SCC(pi) == SCC(¬pi).
Requires `recompute!` first (called automatically).
"""
function has_conflict(net::ImpGraph)
    recompute!(net; compute_reach=false)
    for i in 1:net.n
        if net.comp_id[i] == net.comp_id[i + net.n]
            return true
        end
    end
    return false
end

"""
    conflicts(net) -> Vector{Int}

Return variable indices i where SCC(pi) == SCC(¬pi).
"""
function conflicts(net::ImpGraph)
    recompute!(net; compute_reach=false)
    bad = Int[]
    for i in 1:net.n
        if net.comp_id[i] == net.comp_id[i + net.n]
            push!(bad, i)
        end
    end
    return bad
end

"""
    implies(net, x, y) -> Bool

Check whether x ⇒ y holds in the SCC DAG (i.e., path SCC(x) → SCC(y)).
Uses reachability cache if available; otherwise falls back to `has_path`.
"""
function implies(net::ImpGraph, x::Int, y::Int)
    recompute!(net; compute_reach=true)
    cx = net.comp_id[x]
    cy = net.comp_id[y]
    if !isempty(net.reach)
        return net.reach[cx][cy]
    else
        return has_path(net.cg, cx, cy)
    end
end

"""
    forced_true_literals(net) -> Vector{Int}

Implements your slide fixing rule:
If for any literal v we have SCC(¬v) → SCC(v), then all literals in SCC(v) can be fixed true.

Returns a list of literal node ids that are forced true.
(You can translate to variable fixings depending on your conventions.)
"""
function forced_true_literals(net::ImpGraph)
    recompute!(net; compute_reach=true)
    forced = BitVector(falses(num_literals(net)))

    for v in 1:num_literals(net)
        cv = net.comp_id[v]
        cnot = net.comp_id[neg(net, v)]
        if !isempty(net.reach) && net.reach[cnot][cv]
            # all nodes in SCC(cv) are forced true
            for w in net.comps[cv]
                forced[w] = true
            end
        end
    end

    return findall(forced)
end

"""
    forced_variable_assignments(net) -> Dict{Int,Bool}

Return forced assignments for variables:
- if literal pi forced true => i => true
- if literal ¬pi forced true => i => false
If both appear, you already have a conflict.
"""
function forced_variable_assignments(net::ImpGraph)
    lits = forced_true_literals(net)
    asg = Dict{Int,Bool}()
    for v in lits
        (i, isneg) = var_of(net, v)
        val = !isneg
        if haskey(asg, i) && asg[i] != val
            # conflicting forced assignment; keep both info by throwing
            throw(ErrorException("Conflict in forced assignments for p$i"))
        end
        asg[i] = val
    end
    return asg
end

# ---------------------------
# Pretty helpers
# ---------------------------

"""
    literal_string(net, v)

Human-readable: "p7" or "¬p7".
"""
function literal_string(net::ImpGraph, v::Int)
    (i, isneg) = var_of(net, v)
    return isneg ? "¬p$i" : "p$i" 
end