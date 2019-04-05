struct SimpleNet
    S::Array{Float64, 2}
    b::Array{Float64,1}
    lb::Array{Float64,1}
    ub::Array{Float64,1}
    idxflux::Array{Int,1}
    idxmeta::Array{Int,1}
end

function reduceiterative(X::COBRA.LPproblem; tiny::Real=0.0)
    RX = deepcopy(X)
    retvarval = Dict{String,Float64}()
    while true
        varval, Xnew = findminmax(RX,tiny=tiny)
        nfixed = length(varval)
        println("number of fixed variables = ", length(varval))
        if nfixed == 0 
            return retvarval, RX
        end
        for i=1:nfixed
           # println(retvarval," nfixed = $nfixed ", varval[i])
            fluxname = Xnew.rxns[varval[i][1]]
            haskey(retvarval, fluxname) && warn("key $fluxname already present")
            retvarval[fluxname] = varval[i][2]
        end
        scraX, RX = reduceproblem(Xnew, varval)
        println("size of RX = ", size(RX.S))
#        return retvarval,RX
    end
end

function findminmax(X::COBRA.LPproblem;
                    solverName::Symbol=:Gurobi,
                    solParams=[],
                    optPercentage::Float64=100.0,
                    tiny::Real=1e-10)

    solver = COBRA.changeCobraSolver(solverName, solParams)    
    Xnew = deepcopy(X)           
    minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax = COBRA.distributedFBA(X, solver, nWorkers=1, optPercentage=optPercentage)       
    varval = Array{Tuple{Int,Float64},1}()
    for i in eachindex(maxFlux)
        deltaflux = maxFlux[i] - minFlux[i] 
        if -tiny <= deltaflux <= tiny 
            if abs(maxFlux[i]) < tiny
                maxFlux[i] = 0.0
                minFlux[i] = 0.0
            else
                maxFlux[i] = 0.5*(minFlux[i] + maxFlux[i])
                minFlux[i] = maxFlux[i]
            end
            push!(varval,(i, minFlux[i]))
        end        
    end    
    return varval, COBRA.LPproblem(X.S,X.b,X.c,minFlux,maxFlux,X.osense,X.csense,X.rxns,X.mets)
end

function reduceproblem(X::COBRA.LPproblem,fixedvar::Array{Tuple{Int64,Float64},1})
    res = _reduceproblem(X.S,X.b,X.lb,X.ub,fixedvar)
    @extract res : S b lb ub idxflux idxmeta
    res, COBRA.LPproblem(sparse(S), b, X.c[idxflux], lb, ub, X.osense, X.csense[idxmeta], X.rxns[idxflux], X.mets[idxmeta])
end

function _reduceproblem(S, b, lb, ub, fixedvar::Array{Tuple{Int64,Float64},1})    
    M,N0 = size(S)
    fS0 = full(S)

    idxfix = collect(map(x->x[1],fixedvar))
    valfix = collect(map(x->x[2],fixedvar))

    idxnonfix = setdiff(1:N0,idxfix) # indices of non-fixed variables 

    N = length(idxnonfix) # number variables reduced problem
    nS=zeros(Float64, M, N)
    nlb=zeros(Float64, N)
    nub=zeros(Float64, N)
    nb=zeros(Float64, M) 

    deltab = zeros(M)

    for j=1:length(idxfix)
        for i=1:M
            deltab[i] += fS0[i,idxfix[j]]*valfix[j]
        end
    end

    for i=1:M
        nb[i] = b[i]-deltab[i]
    end

    
    for i=1:N
        nlb[i] = lb[idxnonfix[i]]
        nub[i] = ub[idxnonfix[i]]
        for j=1:M
            nS[j,i] = fS0[j,idxnonfix[i]]
        end
    end

    idx = idxlicols(nS')
    nb1 = zeros(length(idx))
    for i=1:length(idx)
        nb1[i] = nb[idx[i]]
    end
    SimpleNet(nS[idx,:], nb1, nlb, nub, idxnonfix, idx)
end
