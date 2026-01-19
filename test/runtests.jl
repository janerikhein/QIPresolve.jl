using QIPresolve
using Test

@testset "QIPresolve.jl" begin
    include("xor_elim_tests.jl")
    include("quadexpr_tests.jl")
end
