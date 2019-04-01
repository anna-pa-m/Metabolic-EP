# utilities that wait for COBRA.jl to update to 1.0
# not included from MetabolicEP.jl

# methods originally in ep.jl

metabolicEP(X::COBRA.LPproblem; args...) = metabolicEP(X.S, X.b, X.lb, X.ub; args...) # convinence for runnning directly on Lpproblems

# methods from here were formerly in utils.jl

function reduceModel(X::COBRA.LPproblem; solverName::Symbol=:Gurobi,solParams=[],optPercentage::Float64=100.0, tiny=1e-10)

    solver = COBRA.changeCobraSolver(solverName, solParams)
    minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax = COBRA.distributedFBA(X, solver, nWorkers=1, optPercentage=optPercentage)    
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
        end
    end
    return COBRA.LPproblem(X.S,X.b,X.c,minFlux,maxFlux,X.osense,X.csense,X.rxns,X.mets)
end



function reduceModel!(X::COBRA.LPproblem; solverName::Symbol=:Gurobi,solParams=[],optPercentage::Float64=100.0, tiny=1e-10)
        
    solver = COBRA.changeCobraSolver(solverName, solParams)
    minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax = COBRA.distributedFBA(X, solver, nWorkers=1, optPercentage=optPercentage)    

    for i in eachindex(X.lb)
        deltaflux = maxFlux[i] - minFlux[i] 
        if deltaflux > tiny 
            X.lb[i] = minFlux[i]
            X.ub[i] = maxFlux[i]
        elseif -tiny <= deltaflux <= tiny 
            if abs(X.lb[i]) < tiny
                X.lb[i] = 0.0
                X.ub[i] = 0.0
            else
                X.lb[i] = 0.5*(minFlux[i] + maxFlux[i])
                X.ub[i] = X.lb[i]
            end
        else
            warn("lb[$i] = ", lb[i], " ub[$i] = ",ub[i])
        end
    end
    return nothing
end

function simplereducestandardform(X)    
    RX = reduceModel(X)  
    @extract RX : S b c lb ub csense osense rxns mets
    M,N = size(S)
    idxfixed = find(lb .== ub .* lb .== 0) 
    idxrow = setdiff(1:M,idxfixed)
    idxcol = setdiff(1:N,idxfixed)
    COBRA.LPproblem(S[idxrow,idxcol],b[idxrow],c[idxcol],lb[idxcol],ub[idxcol],osense, csense[idxrow], rxns[idxcol], mets[idxrow])       
end

function standardform(X::COBRA.LPproblem)
    idxrow, idxcol, res, bnew = echelonize(full(X.S), X.b)
    (length(idxrow),length(idxcol)) == size(res) || BoundsError()    
    COBRA.LPproblem(sparse(res),bnew,X.c[idxcol],X.lb[idxcol],X.ub[idxcol],X.osense, X.csense[idxrow], X.rxns[idxcol], X.mets[idxrow])    
end

isstandardform(X::COBRA.LPproblem) = isstandardform(X.S)

