using Test

@testset "QIPresolve.jl" begin
    include("quadexpr_tests.jl")
    include("constraint_tests.jl")
    include("interaction_graph_tests.jl")
    include("tree_decomposition_tests.jl")
    include("dp_propagation_tests.jl")
    include("model_tests.jl")
    include("model_stats_tests.jl")
    include("model_io_tests.jl")
    include("instance_generation_tests.jl")
    include("propagation_tests.jl")
    include("xor_constraint_tests.jl")
    include("xor_model_tests.jl")
    include("presolve_parity_tests.jl")
    include("parity_stats_tests.jl")
    include("residue_stats_tests.jl")
    include("presolve_tests.jl")
end
