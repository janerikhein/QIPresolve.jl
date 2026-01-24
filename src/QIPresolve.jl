module QIPresolve

include("core/PresolvingCore.jl")
include("graph_embedding/GraphEmbedding.jl")
include("model_io/ModelIO.jl")
include("config.jl")


import .PresolvingCore
import .GraphEmbedding
import .PresolvingIO
import .PresolveConfig


end
