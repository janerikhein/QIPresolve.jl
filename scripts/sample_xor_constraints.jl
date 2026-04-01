#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

import QIPresolve.PresolvingCore as PC

function symmetric_bitmatrix(n::Int, edges::Vector{Tuple{Int, Int}})
    mat = falses(n, n)
    for (i, j) in edges
        mat[i, j] = true
        mat[j, i] = true
    end
    return mat
end

function linear_terms(con::PC.XorConstraint)
    return findall(con.par)
end

function conjunction_terms(con::PC.XorConstraint)
    con.conj === nothing && return Tuple{Int, Int}[]
    pairs = Tuple{Int, Int}[]
    for i in 1:size(con.conj, 1)
        for j in (i + 1):size(con.conj, 2)
            con.conj[i, j] && push!(pairs, (i, j))
        end
    end
    return pairs
end

function constraint_holds(con::PC.XorConstraint, x::AbstractVector{Bool})
    lhs = false

    for i in eachindex(con.par)
        con.par[i] || continue
        lhs = xor(lhs, x[i])
    end

    if con.conj !== nothing
        for i in 1:length(x)
            xi = x[i]
            for j in (i + 1):length(x)
                con.conj[i, j] || continue
                lhs = xor(lhs, xi & x[j])
            end
        end
    end

    return lhs == con.rhs
end

function show_constraint(title::String, con::PC.XorConstraint)
    println(title)
    println("  linear terms:       ", linear_terms(con))
    println("  conjunction terms:  ", conjunction_terms(con))
    println("  rhs:                ", con.rhs)
    println(
        "  meta:               ",
        (
            nnz_par = con.meta.nnz_par,
            nnz_conj = con.meta.nnz_conj,
            is_pure_xor = con.meta.is_pure_xor,
            requires_update = con.meta.requires_update,
            requires_prop = con.meta.requires_prop,
        ),
    )
    println()
    return nothing
end

println("XorConstraint demo")
println("==================")
println()

con1 = PC.XorConstraint(
    BitVector([1, 0, 1, 0]),
    symmetric_bitmatrix(4, [(1, 2), (2, 4)]),
    true,
)
show_constraint("1. Mixed XOR/XOR-AND constraint", con1)

con2 = PC.XorConstraint(
    BitVector([1, 1, 0, 0]),
    symmetric_bitmatrix(4, [(2, 4)]),
    false,
)
PC.xor_con!(con1, con2)
show_constraint("2. After xor_con!(con1, con2)", con1)
PC.update!(con1)
show_constraint("3. After update!(con1)", con1)

bipartite = PC.XorConstraint(
    falses(4),
    symmetric_bitmatrix(4, [(1, 3), (1, 4), (2, 3), (2, 4)]),
    true,
)
show_constraint("4. Full bipartite conjunction row", bipartite)
split = PC.split_bipartite(bipartite)
if split === nothing
    println("split_bipartite returned nothing")
    println()
else
    left, right = split
    show_constraint("5. split_bipartite left row", left)
    show_constraint("6. split_bipartite right row", right)
end

original_single = PC.XorConstraint(
    BitVector([1, 0, 0, 1]),
    symmetric_bitmatrix(4, [(1, 2), (1, 3), (2, 4)]),
    false,
)
single_sub = deepcopy(original_single)
PC.substitute_var!(single_sub, 1, 2, true)
PC.update!(single_sub)
show_constraint("7. After substitute_var!(..., 1, 2, true)", single_sub)

x_single = BitVector([1, 0, 1, 0])
x_single[1] = xor(x_single[2], true)
println("  sample assignment for single substitution: ", collect(x_single))
println("  original holds:    ", constraint_holds(original_single, x_single))
println("  transformed holds: ", constraint_holds(single_sub, x_single))
println()

original_mask = PC.XorConstraint(
    BitVector([1, 0, 1, 0, 0]),
    symmetric_bitmatrix(5, [(1, 2), (1, 4), (1, 5)]),
    false,
)
mask_sub = deepcopy(original_mask)
mask = BitVector([0, 1, 0, 1, 0])
PC.substitute_var!(mask_sub, 1, mask, true)
PC.update!(mask_sub)
show_constraint("8. After mask-based substitution", mask_sub)

x_mask = BitVector([0, 1, 1, 0, 1])
x_mask[1] = xor(true, xor(x_mask[2], x_mask[4]))
println("  sample assignment for mask substitution: ", collect(x_mask))
println("  original holds:    ", constraint_holds(original_mask, x_mask))
println("  transformed holds: ", constraint_holds(mask_sub, x_mask))
