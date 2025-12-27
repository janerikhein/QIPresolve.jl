using Base.Threads

using BenchmarkTools


function xor_vec_simd(X)
    @inbounds @simd for i=2:length(X)
        X[i] .⊻= X[1]
    end
end


function xor_vec(X)
    @inbounds for i=2:length(X)
        X[i] .⊻= X[1]
    end
end

function xor_mat_simd(X)
    @inbounds @simd for j=2:size(X,2)
        @inbounds @simd for i=1:size(X,1) 
            X[i, j] ⊻= X[i, 1]
        end
    end
end


#@benchmark xor_vec(X) setup=(X=[BitVector(rand(Bool, 64*10^2)) for i=1:10^4])
#@benchmark xor_mat_simd(X) setup=(X=BitMatrix(rand(Bool, 64*10^2, 10^4)))


#@benchmark xor_vec_th(X) setup=(X=[BitVector(rand(Bool, 64*10^2)) for i=1:10^4])

'''
mutable struct DenseXORModel
    rows::Vector{BitVector}
    perm::Vector{Int}

    function DenseXORModel(mat::Vector{BitVector})
        perm = collect(1:size(mat, 1))
        new(mat, perm)
    end

end

function DenseXORModel(nrows::Int, ncols::Int)
    rows = Vector{BitVector}(undef, nrows)

    for i in 1:nrows
        rows[i] = BitVector(undef, ncols)
    end

    DenseXORModel(rows)
end
'''