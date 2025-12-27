module GaussJordanGF2

export FirstOnePivot, gauss_jordan_gf2!

"""
    PivotRule

Abstract policy for selecting pivots in Gauss–Jordan elimination over GF(2).
"""
abstract type PivotRule end

"""
    FirstOnePivot <: PivotRule

Pivoting rule that selects the next non-zero row in logical order, and within that row, the leftmost non-zero column.

"""
struct FirstOnePivot <: PivotRule end

"""
    PivotState

Abstract container type for caching data for pivoting rules.
"""
abstract type PivotState end

"""
    FirstOneState

Pivot state for 'FirstOnePivot' rule. Stores the index 'rstart' of the next row to consider
"""
mutable struct FirstOneState <: PivotState
    rstart::Int
end

initialize_pivot_state(::FirstOnePivot)::FirstOneState = FirstOneState(1)

function get_pivot!(state::FirstOneState, rows::Vector{BitVector})
    
    # ignore RHS entry in rows
    nvars = length(rows[1]) - 1
    m = length(rows)

    @inbounds for r in state.rstart:m
        row = rows[r]

        # find first non-zero coefficient in the row
        p = findnext(row, 1)

        if p !== nothing && p <= nvars
            state.rstart = r + 1
            return (r, p)
        end
        # else: row is zero on coefficients → skip it
    end

    return (0, 0)
end


@inline xor_rows!(rows::Vector{BitVector}, source::Int, target::Int) = (rows[target] .⊻= rows[source])

function gauss_jordan_gf2!(rows::Vector{BitVector}, pivot_rule::T) where {T<:PivotRule}

    pivot_state = initialize_pivot_state(pivot_rule)

    while true

        #identify pivot element
        pivot_row_idx, pivot_col_idx = get_pivot!(pivot_state, rows)
        
        if pivot_row_idx == 0 && pivot_col_idx == 0
            break
        end
        
        # eliminate all entries in pivot column in all other rows
        @inbounds for row_idx in eachindex(rows)
            if row_idx === pivot_row_idx
                continue
            end
            if rows[row_idx][pivot_col_idx]
                xor_rows!(rows, pivot_row_idx, row_idx)
            end
        end
        
    end
end


end # module




