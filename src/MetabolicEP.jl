module MetabolicEP
using COBRA

@static VERSION > v"0.6" && using SpecialFunctions

export metabolicEP, ReadMatrix, MetNet, EPFields, EPout

include("ep.jl")
include("utils.jl")

end # end module

