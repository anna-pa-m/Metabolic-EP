using JuMP, Gurobi, SCS, GLPK,  ExtractMacro, SparseArrays
mutable struct ReducedNet
    S::Matrix{Float64}
    b::Vector{Float64}
    biter::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    idxflux::Vector{Int}
    idxmeta::Vector{Int}
end

ReducedNet(S,b,lb,ub)=ReducedNet(S,b,b,lb,ub,[1:size(S,2);],MetabolicEP.idxlicols(S'))
ReducedNet(S::AbstractSparseArray,b,lb,ub)=ReducedNet(Matrix(S),b,lb,ub) 

function reduceiterative!(X::ReducedNet; tiny::Real=0.0)
    @extract X S biter
    nfixed = 0
    RX = deepcopy(X)
    varval = Dict{Int,Float64}()
    while true
        jumpminmax!(RX)
        fixfluxes!(RX,varval,tiny)
        # Xnew = ReducedNet(S,RX.b,minFlux,maxFlux,RX.idxflux)         
        if length(varval) == nfixed 
            return varval, RX
        end
        reduceproblem!(RX, varval)
        nfixed = length(varval)
        println("fixed $nfixed variables. Free variable = ", length(RX.idxflux))
    end
end


function jumpminmax!(X::ReducedNet;
                    solverName::Symbol=:Gurobi,
                    verbose=0,
                    solParams=[],
                    optPercentage::Float64=100.0,
                    tiny::Real=1e-10)

    @extract X: S biter lb ub idxflux idxmeta 

    vS = view(S,idxmeta,idxflux)
    vlb = view(lb,idxflux)
    vub = view(ub,idxflux)
    vbiter = view(biter,idxmeta)
    
    N = length(idxflux)
    model = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=verbose))
    @variable(model, x[i=1:N], lower_bound=vlb[i],upper_bound=vub[i])
    @constraint(model, vS*x .== vbiter) 
    for l in 1:N
        @objective(model, Min , x[l]) # minimizer 
        optimize!(model)
        termination_status(model) == MOI.OPTIMAL || (@warn "problem with min of $l variable")
        vlb[l] = objective_value(model)
        @objective(model, Max , x[l]) # maximizer
        optimize!(model)
        termination_status(model) == MOI.OPTIMAL || (@warn "problem with max of $l variable")
        vub[l] = objective_value(model)               
    end
end

function fixfluxes!(X::ReducedNet, varval::Dict{Int,Float64},  tiny::Float64)
  
    @extract X: lb ub idxflux
    for i in idxflux
        deltaflux = ub[i] - lb[i] 
        if -tiny <= deltaflux <= tiny 
            if max(abs(ub[i]),abs(lb[i])) < tiny
                lb[i] = 0.0
                ub[i] = 0.0
            else
                ub[i] = (lb[i] + ub[i])/2
                lb[i] = ub[i]
            end
            if haskey(varval,i) && varval[i] != lb[i]
                @warn("varval[$i] was $(varval[i]) now I want to set it to  $(lb[i])")
            else
                varval[i] = lb[i]
            end
        end
    end   
end


function reduceproblem!(X::ReducedNet, varval::Dict{Int,Float64})

    @extract X: S b biter lb ub idxflux idxmeta
    M,N0 = size(S)
    idxfix = collect(keys(varval))
    valfix = collect(values(varval))
    idxnonfix = setdiff(1:N0,idxfix) # indices of non-fixed variables 
    X.idxflux=idxnonfix
    N = length(idxnonfix) # number variables reduced problem
    # println("N=$N ", length(idxflux))           
    # nS=zeros(Float64, M, N)
    # nlb=zeros(Float64, N)
    # nub=zeros(Float64, N)
    # nb=zeros(Float64, M) 
    deltab = zero(b)
    for j ∈ eachindex(idxfix)
        for i ∈ idxmeta
            deltab[i] += S[i,idxfix[j]]*valfix[j]
        end
    end
    
    for i ∈ idxmeta
        biter[i] -= deltab[i]
    end

    vS = view(S,:,idxnonfix)
    X.idxmeta = MetabolicEP.idxlicols(vS')
    nothing
end

# function reduceproblem(X::COBRA.LPproblem,fixedvar::Array{Tuple{Int64,Float64},1})
#     res = _reduceproblem(X.S,X.b,X.lb,X.ub,fixedvar)
#     @extract res : S b lb ub idxflux idxmeta
# end

# Function findminmax(X::COBRA.LPproblem;
#                     solverName::Symbol=:Gurobi,
#                     solParams=[],
#                     optPercentage::Float64=100.0,
#                     tiny::Real=1e-10)

#     solver = COBRA.changeCobraSolver(solverName, solParams)    
#     Xnew = deepcopy(X)       
#     minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax = COBRA.distributedFBA(X, solver, nWorkers=1, optPercentage=optPercentage)       
#     varval = Array{Tuple{Int,Float64},1}()
#     for i in eachindex(maxFlux)
#         deltaflux = maxFlux[i] - minFlux[i] 
#         if -tiny <= deltaflux <= tiny 
#             if abs(maxFlux[i]) < tiny
#                 maxFlux[i] = 0.0
#                 minFlux[i] = 0.0
#             else
#                 maxFlux[i] = 0.5*(minFlux[i] + maxFlux[i])
#                 minFlux[i] = maxFlux[i]
#             end
#             push!(varval,(i, minFlux[i]))
#         end        
#     end    
#     return varval, COBRA.LPproblem(X.S,X.b,X.c,minFlux,maxFlux,X.osense,X.csense,X.rxns,X.mets)
# end

# function reduceproblem(X::COBRA.LPproblem,fixedvar::Array{Tuple{Int64,Float64},1})
#     res = _reduceproblem(X.S,X.b,X.lb,X.ub,fixedvar)
#     @extract res : S b lb ub idxflux idxmeta
#     res, COBRA.LPproblem(sparse(S), b, X.c[idxflux], lb, ub, X.osense, X.csense[idxmeta], X.rxns[idxflux], X.mets[idxmeta])
# end

