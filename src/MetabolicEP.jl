module MetabolicEP
#using COBRA, MAT, ExtractMacro
using MAT, ExtractMacro, SpecialFunctions, SparseArrays
using SparseArrays: SparseMatrixCSC
using Printf: @printf
using Clp, MathProgBase

include("HitAndRun.jl")
using .HitAndRun # internal module

VERSION >= v"0.6.0-rc1" && using SpecialFunctions

export metabolicEP, ReadMatrix, MetNet, EPFields, EPout #, standardform, reduceModel, reduceiterative
export HitAndRun: hrsample
include("types.jl")
include("ep.jl")
include("utils.jl")

# we wait for COBRA to update ...

# include("reduceiterative.jl")

end # end module
