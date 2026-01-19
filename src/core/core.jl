module PresolvingCore

include("qp_model.jl")
include(".basic_reductions/basic_reductions.jl")
include(".parity_reductions/parity_reductions.jl")
include(".residue_reductions/residue_reductions.jl")

using .BasicReductions
using .ParityReductions
using .ResidueReductions

end # module