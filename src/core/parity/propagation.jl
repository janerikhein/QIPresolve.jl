using DataStructures: Deque

const VarId = Int

@enum LitLabel TRUE FALSE UNDEF

struct VarLit
    vid::VarId
    neg::Bool
end

negated(lit::VarLit) = VarLit(lit.vid, !lit.neg)

struct SCCNode
    lits::BitVector
    lits_neg::BitVector
end

mutable struct PropagationManager
    nvars::Int
    nreps::Int

    # Literal index mappings for UnionFind structure
    lit_to_pos::Dict{VarLit, Int}
    pos_to_lit::Vector{VarLit}

    # UnionFind tree structure (parent_pos[i] == i, iff i is representant) 
    parent_pos::Vector{Int}

    # Index mapping for implication graph storing only representant literals
    rep_pos_to_scc_pos::Dict{Int, Int}
    scc_pos_to_rep_pos::Vector{Int}

    # Labels and deques for querying pending fixations and substitutions 
    var_labels::Dict{VarId, Bool}
    

    # Implication graph structure over SCC representant positions
    graph::DiGraph 
end


function repr!(manager::PropagationManager, pos::Int)
    root = pos
    while prop.parent_pos[root] != root
        root = prop.parent_pos[root]
    end

    # path compression
    while pos != root
        p = prop.parent_pos[pos]
        prop.parent_pos[pos] = root
        pos = p
    end

    return root
end

function union_reps!(manager::PropagationManager, rep_positions::Vector{Int})
    new_rep = min(rep_positions)

    for pos in reps
        @assert manager.parent_pos[pos] == pos
        manager.parent_pos[pos] = new_rep
    end

    return new_rep
end


function add_edge!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)
    rep1 = repr!(manager.lit_to_pos[lit1])
    rep2 = repr!(manager.lit_to_pos[lit2])
    added = false
    
    if rep1 != rep2
        scc_pos1 = manager.rep_pos_to_scc_pos[rep1]
        scc_pos2 = manager.rep_pos_to_scc_pos[rep2]
        added = add_edge!(manager.graph, scc_pos1, scc_pos2)
    end

    return added
end


function add_implication!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)
    add_edge!(manager, lit1, lit2)
    add_edge!(manager, negated(lit2), negated(lit1))
end


function add_equivalence!(manager::PropagationManager, lit1::VarLit, lit2::VarLit)
    add_implication!(manager, lit1, lit2)
    add_implication!(manager, lit2, lit1)
end


function fix_var!(manager::PropagationManager, vid::VarId, val::Bool)
    lit = VarLit(vid, val)
    add_edge!(manager, negated(lit), lit)
end


function update_sccs!(manager::PropagationManager)

    # get all current representants grouped into new sccs (transformed to lit positions)
    comps = strongly_connected_components(manager.graph) .|> (arr -> manager.scc_pos_to_rep_pos[arr])
    ncomps = length(comps)

    # union comps and update scc mappings
    new_rep_to_idx = Dict{Int, Int}
    sizehint!(new_rep_to_idx, ncomps)
    new_idx_to_rep = Vector{Int}(undef, ncomps)

    for (i, comp) in enumerate(comps)
        new_rep = union!(manager, comp)
        new_rep_to_idx[new_rep] = i
        new_idx_to_rep[i] = new_rep
    end
    
    # build new condensation graph
    new_graph = DiGraph(ncomps)
    for edge in edges(graph)
        src_rep = repr!(manager.scc_pos_to_rep_pos[src(edge)])
        dst_rep = repr!(manager.scc_pos_to_rep_pos[dst(edge)])
        if src_rep != dst_rep
            src_idx = new_rep_to_idx[src_rep]
            dst_idx = new_rep_to_idx[dst_idx]
            add_edge!(new_graph, src_idx, dst_idx)
        end
    end  

    return manager
end

