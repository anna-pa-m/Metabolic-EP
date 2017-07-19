module MetabolicEP
using COBRA, MAT, ExtractMacro

VERSION >= v"0.6.0-rc1" && using SpecialFunctions

export metabolicEP, ReadMatrix, MetNet, EPFields, EPout

include("types.jl")
include("ep.jl")
include("utils.jl")

end # end module

