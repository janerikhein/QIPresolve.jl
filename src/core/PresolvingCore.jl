module PresolvingCore

export
    # core model functionality
    QPModel,
    register_con!,
    register_var_info!,
    register_obj!
    

include("qp_model.jl")
include("model_io.jl")
include("basic_reductions/basic_reductions.jl")
include("parity_reductions/parity_reductions.jl")
include("residue_reductions/residue_reductions.jl")

end # module
