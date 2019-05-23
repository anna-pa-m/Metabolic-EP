using MetabolicEP,Test

tests = ["testep"]

printstyled("Running tests:\n", color=:blue)

using Random
Random.seed!(345679)

res = map(tests) do t
    @eval module $(Symbol("Test_",t))
    using Random
    Random.seed!(345679)
    include($t * ".jl")
    end
    return
end
# print method ambiguities
println("Potentially stale exports: ")
ambiguities=Test.detect_ambiguities(PottsGauge)
isempty(ambiguities) || display(ambiguities)
println()
