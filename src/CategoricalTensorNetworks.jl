module CategoricalTensorNetworks

using Reexport

include("TensorNetworkAlgebras.jl")
include("ScheduleUWDs.jl")

@reexport using .TensorNetworkAlgebras
@reexport using .ScheduleUWDs

end
