using Test

@testset "QIPresolve.jl" begin
    include("quadexpr_tests.jl")
    include("constraint_tests.jl")
    include("model_tests.jl")
end
