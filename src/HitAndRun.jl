module HitAndRun
using LinearAlgebra
using SparseArrays, JuMP, GLPK
export hrsample, warmup
include("utils.jl")

function hrsample(S, b, lb, ub; niter=10^6, nsamples = 10000)
    x = warmup(S,b,lb,ub);
    m,n = size(S)
    # preprocessing: find base
    base = nullspace(Matrix(S))
    k = size(base,2)
    v = zeros(n)
    h = zeros(nsamples,n)
    hpos = 1
    for it=1:niter
        # pick a random direction
        v .= base * randn(k)
        dx = v/norm(v)
        # compute intersection
        l = maximum(min((lb[i]-x[i])/dx[i], (ub[i]-x[i])/dx[i]) for i=1:n if dx[i] != 0)
        u = minimum(max((lb[i]-x[i])/dx[i], (ub[i]-x[i])/dx[i]) for i=1:n if dx[i] != 0)
        # find a random point in the intersection
        t = l+(u-l)*rand()
        x .+= t * dx
        if mod(it, floor(niter/nsamples)) == 0
            h[hpos,:] .= x[:]
            hpos += 1
        end
    end
    h
end

function warmup(S, b, lb, ub)
    n = length(lb)
    x0 = zeros(n)
    ei = zeros(n)
    for i=1:n
        ei[i] = -1.0
        sol=linprog(ei, S, fill('=', length(b)), b, lb, ub, GLPK.Optimizer)
        x0 .+= sol.sol / 2n
        ei[i] = +1.0
        sol=linprog(ei, S, fill('=', length(b)), b, lb, ub, GLPK.Optimizer)
        ei[i] = 0.0
        x0 .+= sol.sol / 2n
    end
    x0
end
end #end module HitAndRun