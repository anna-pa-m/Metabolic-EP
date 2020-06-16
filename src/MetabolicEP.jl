module MetabolicEP
#using COBRA, MAT
using MAT, ExtractMacro, SpecialFunctions, SparseArrays, Printf
using SparseArrays: SparseMatrixCSC
using Printf: @printf
using Clp, MathProgBase

include("HitAndRun.jl")
using .HitAndRun # internal module

export metabolicEP, ReadMatrix, MetNet, EPFields, EPout #, standardform, reduceModel, reduceiterative
export hrsample
include("types.jl")
include("ep.jl")
include("utils.jl")

# we wait for COBRA to update ...

# include("reduceiterative.jl")

end # end module
