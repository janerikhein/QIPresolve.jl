import Graphs
using Graphs: SimpleGraph

"""
    InteractionComponent

Abstract supertype for separable components induced by an [`InteractionGraph`](@ref).
"""
abstract type InteractionComponent end

"""
    LinSingleton(var_id, lin_coeff)

Single-variable linear residual component.
"""
struct LinSingleton <: InteractionComponent
    var_id::VarId
    lin_coeff::Int
end

"""
    QuadSingleton(var_id, lin_coeff, quad_coeff)

Single-variable quadratic residual component with a nonzero diagonal term.
"""
struct QuadSingleton <: InteractionComponent
    var_id::VarId
    lin_coeff::Int
    quad_coeff::Int
end

"""
    NonSingleton(graph, var_id_to_pos, pos_to_var_id, lin_coeffs, quad_coeffs)

Multi-variable residual component induced by one connected component of an
interaction graph.
"""
struct NonSingleton <: InteractionComponent
    graph::SimpleGraph{Int}
    var_id_to_pos::Dict{VarId, Int}
    pos_to_var_id::Vector{VarId}
    lin_coeffs::Dict{VarId, Int}
    quad_coeffs::Dict{Tuple{VarId, VarId}, Int}
end

"""
    InteractionGraph(con, m)

Build the bilinear interaction graph induced by one quadratic constraint modulo
`m`.

Vertices correspond to active variables in `con.qe`, ordered by variable id.
An undirected edge is present for each off-diagonal quadratic coefficient that
does not reduce to zero modulo `m`.
"""
struct InteractionGraph
    graph::SimpleGraph{Int}
    var_id_to_pos::Dict{VarId, Int}
    pos_to_var_id::Vector{VarId}
    modulus::Int
    lin_coeffs::Dict{VarId, Int}
    quad_coeffs::Dict{Tuple{VarId, VarId}, Int}
end

function _integer_residue(coeff::Float64, modulus::Int, description::String)
    isinteger(coeff) || throw(ArgumentError("$description must be integer-valued, got $coeff"))
    return mod(trunc(Int, coeff), modulus)
end

function InteractionGraph(con::Constraint, m::Integer)
    m > 0 || throw(ArgumentError("modulus must be positive, got $m"))
    modulus = Int(m)

    pos_to_var_id = sort!(collect(vars(con.qe)))
    var_id_to_pos = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var_id))
    graph = SimpleGraph{Int}(length(pos_to_var_id))
    lin_coeffs = Dict{VarId, Int}()
    quad_coeffs = Dict{Tuple{VarId, VarId}, Int}()

    for i in 1:length(pos_to_var_id)
        vid_i = pos_to_var_id[i]
        lin_residue = _integer_residue(
            get_lin_coeff(con.qe, vid_i),
            modulus,
            "linear coefficient for $vid_i",
        )
        lin_residue == 0 || (lin_coeffs[vid_i] = lin_residue)

        diag_residue = _integer_residue(
            get_quad_coeff(con.qe, vid_i, vid_i),
            modulus,
            "diagonal coefficient for ($vid_i, $vid_i)",
        )
        diag_residue == 0 || (quad_coeffs[(vid_i, vid_i)] = diag_residue)

        for j in (i + 1):length(pos_to_var_id)
            vid_j = pos_to_var_id[j]
            residue = _integer_residue(
                get_quad_coeff(con.qe, vid_i, vid_j),
                modulus,
                "bilinear coefficient for ($vid_i, $vid_j)",
            )
            residue == 0 && continue
            quad_coeffs[(vid_i, vid_j)] = residue
            Graphs.add_edge!(graph, i, j)
        end
    end

    return InteractionGraph(graph, var_id_to_pos, pos_to_var_id, modulus, lin_coeffs, quad_coeffs)
end

function _component_subgraph(ig::InteractionGraph, pos_to_var_id::Vector{VarId})
    graph = SimpleGraph{Int}(length(pos_to_var_id))

    for i in 1:length(pos_to_var_id)
        src_vid = pos_to_var_id[i]
        src_pos = ig.var_id_to_pos[src_vid]
        for j in (i + 1):length(pos_to_var_id)
            dst_vid = pos_to_var_id[j]
            dst_pos = ig.var_id_to_pos[dst_vid]
            Graphs.has_edge(ig.graph, src_pos, dst_pos) || continue
            Graphs.add_edge!(graph, i, j)
        end
    end

    return graph
end

function _non_singleton(ig::InteractionGraph, component::Vector{Int})
    pos_to_var_id = sort!(ig.pos_to_var_id[component])
    var_id_to_pos = Dict(vid => pos for (pos, vid) in enumerate(pos_to_var_id))
    graph = _component_subgraph(ig, pos_to_var_id)
    component_vars = Set(pos_to_var_id)

    lin_coeffs = Dict{VarId, Int}()
    for vid in pos_to_var_id
        residue = get(ig.lin_coeffs, vid, 0)
        residue == 0 || (lin_coeffs[vid] = residue)
    end

    quad_coeffs = Dict{Tuple{VarId, VarId}, Int}()
    for ((vid_i, vid_j), residue) in ig.quad_coeffs
        vid_i in component_vars || continue
        vid_j in component_vars || continue
        quad_coeffs[(vid_i, vid_j)] = residue
    end

    return NonSingleton(graph, var_id_to_pos, pos_to_var_id, lin_coeffs, quad_coeffs)
end

"""
    decompose(ig) -> Vector{InteractionComponent}

Return the separable residual components induced by the connected components of
`ig`.
"""
function decompose(ig::InteractionGraph)::Vector{InteractionComponent}
    components = Graphs.connected_components(ig.graph)
    sort!(components; by = component -> minimum(ig.pos_to_var_id[component]))

    interaction_components = InteractionComponent[]
    for component in components
        if length(component) == 1
            vid = ig.pos_to_var_id[only(component)]
            lin_coeff = get(ig.lin_coeffs, vid, 0)
            quad_coeff = get(ig.quad_coeffs, (vid, vid), 0)
            if quad_coeff != 0
                push!(interaction_components, QuadSingleton(vid, lin_coeff, quad_coeff))
            elseif lin_coeff != 0
                push!(interaction_components, LinSingleton(vid, lin_coeff))
            end
            continue
        end

        push!(interaction_components, _non_singleton(ig, component))
    end

    return interaction_components
end
