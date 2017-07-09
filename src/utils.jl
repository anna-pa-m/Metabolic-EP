using MAT

type MetNet
    N::Int # number of fluxes
    M::Int # number of metabolites
    S::SparseMatrixCSC{Float64,Int} # Stoichiometric matrix M x N sparse
    b::Array{Float64,1} # right hand side of equation  S ν = b 
    c::Array{Float64,1} # reaction index of biomass 
    lb::Array{Float64,1} # fluxes lower bound M elements vector
    ub::Array{Float64,1} # fluxes upper bound M elements vector 
    genes::Array{String,1} # gene names 
    rxnGeneMat::SparseMatrixCSC{Float64,Int} # 
    grRules::Array{String,1} # gene-reaction rule N elements vector of strings (and / or allowed)
    mets::Array{String,1} # metabolites short-name M elements 
    rxns::Array{String,1} # reactions short-name N elements
    metNames::Array{String,1} # metabolites long-names M elements
    metFormulas::Array{String,1} # metabolites formula M elements
    rxnNames::Array{String,1} # reactions long-names N elements
    rev::Array{Bool,1} # reversibility of reactions N elements
    subSystems::Array{String,1} # cellular component of fluxes N elements
end

function ReadMatrix(filename::String)

    X = matread(filename)
    key1=""
    if length(keys(X)) == 1
        key1 = collect(keys(X))[1]
    else
        for i in keys(X)
            if i == "model"
                key1 = i
                break
            end
        end
    end
    S = X[key1]["S"]
    M,N = size(S)
    b = vec(X[key1]["b"])
    c = vec(X[key1]["c"]) 
    lb = vec(X[key1]["lb"])
    ub = vec(X[key1]["ub"])
    if haskey(X[key1],"genes") 
        genes = String[ string(X[key1]["genes"][i]) for i=1:length(X[key1]["genes"])]
    else
        genes = ["NA"]
    end
    if haskey(X[key1], "rxnGeneMat")
        rxnGeneMat = X[key1]["rxnGeneMat"]
    else
        rxnGeneMat = sparse([1],[1],[0.0])
    end
    if haskey(X[key1],"grRules")
        grRules = String[ string(X[key1]["grRules"][i]) for i=1:length(X[key1]["grRules"]) ]
    else
        grRules = ["NA"]
    end
    if haskey(X[key1],"mets")        
        mets = String[ string(X[key1]["mets"][i])  for i=1:length(X[key1]["mets"]) ]        
        if length(unique(mets)) != M 
            mets = String["met$i" for i=1:M]
            warn("not unique list if metabolite names")
        end
    else
        mets = String["met$i" for i=1:M]
    end

    if haskey(X[key1],"rxns")
        rxns = String[ string(X[key1]["rxns"][i])  for i=1:length(X[key1]["rxns"])]
        if length(unique(rxns)) != N 
            rxns = String["rxn$i" for i=1:N]
            warn("not unique list of reaction names")
        end
    else
        rxns = String["rxn$i" for i=1:N]
    end
    if haskey(X[key1],"metNames")
        metNames = String[ string(X[key1]["metNames"][i])  for i=1:length(X[key1]["metNames"]) ]
    else
        metNames = ["NA"]
    end
    if haskey(X[key1],"metFormulas")
        metFormulas = String[ string(X[key1]["metFormulas"][i]) for i=1:length(X[key1]["metFormulas"]) ]
    else
        metFormulas = ["NA"]
    end
    if haskey(X[key1],"rxnNames")
        rxnNames = String[ string(X[key1]["rxnNames"][i])  for i=1:length(X[key1]["rxnNames"]) ]
    else
        rxnNames = ["NA"]
    end
    if haskey(X[key1],"rev")
        rev = Bool[X[key1]["rev"][i] > 0 ? true : false for i=1:length(X[key1]["rev"]) ]
    else
        rev = Vector{Bool}()
    end
    if haskey(X[key1],"subSystems")
        subSystems = String[string(X[key1]["subSystems"][i]) for i=1:length(X[key1]["subSystems"]) ]
    else
        subSystems = ["NA"]
    end
    return MetNet(N, M, S, b, c,  lb, ub, genes, rxnGeneMat, grRules,
                  mets,rxns,metNames,metFormulas, rxnNames,rev,subSystems)

end


include("$(Pkg.dir("COBRA"))/config/solverCfg.jl")

function reduceModel(X::COBRA.LPproblem; solverName::Symbol=:Gurobi,solParams=[],optPercentage::Float64=100.0)

    solver = COBRA.changeCobraSolver(solverName, solParams)
    minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax = COBRA.distributedFBA(X, solver, nWorkers=1, optPercentage=optPercentage)    

    return COBRA.LPproblem(X.S,X.b,X.c,minFlux,maxFlux,X.osense,X.csense,X.rxns,X.mets)

end

function reduceModel!(X::COBRA.LPproblem; solverName::Symbol=:Gurobi,solParams=[],optPercentage::Float64=100.0)

    solver = COBRA.changeCobraSolver(solverName, solParams)
    minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax = COBRA.distributedFBA(X, solver, nWorkers=1, optPercentage=optPercentage)    

    for i in eachindex(X.lb)
        X.lb[i] = minFlux[i]
        X.ub[i] = maxFlux[i]
    end
    return nothing
end
