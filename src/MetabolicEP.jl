module MetabolicEP
using COBRA, MAT, ExtractMacro

VERSION >= v"0.6.0-rc1" && using SpecialFunctions

export metabolicEP, ReadMatrix, MetNet, EPFields, EPout, standardform, reduceModel, reduceiterative

include("types.jl")
include("ep.jl")
include("utils.jl")
include("reduceiterative.jl")

end # end module

