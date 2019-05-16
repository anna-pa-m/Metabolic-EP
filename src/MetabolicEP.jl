module MetabolicEP
#using COBRA, MAT, ExtractMacro
using MAT, ExtractMacro, SpecialFunctions, SparseArrays
using SparseArrays: SparseMatrixCSC
using Printf: @printf
using Clp, MathProgBase


VERSION >= v"0.6.0-rc1" && using SpecialFunctions

export metabolicEP, ReadMatrix, MetNet, EPFields, EPout #, standardform, reduceModel, reduceiterative

include("types.jl")
include("ep.jl")
include("utils.jl")

# we wait for COBRA to update ...

# include("reduceiterative.jl")

end # end module

