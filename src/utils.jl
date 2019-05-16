eye(N) = Matrix(1.0I, N, N)
speye(N) = sparse(eye(N))


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


# function idxlicols(X; tol::Float64=1e-10)
#     sum(abs2,X) == 0 && (return(Array{Int,1}(),Array{Int,2}()))
#     Q,R,E = qr(X,Val{true})  
#     diagr = abs.(diag(R))
#     r = find(diagr .>= tol*diagr[1])[end]
#     idx = sort(E[1:r])
#     return idx
# end

isstandardform(S::SparseMatrixCSC) = sum(abs2,view(S,1:size(S,1),1:size(S,1)) - I) < 1e-12
isstandardform(S::DenseMatrix) = sum(abs2,S[1:size(S,1),1:size(S,1)] -eye(size(S,1))) < 1e-12

function idxlicols(X; tol::Float64=1e-10)
#    sum(abs2,X) == 0 && (return(Array{eltype(X),2}()))
#    res = qr(X,Val(true))
#    return res.p
    sum(abs2,X) == 0 && (return(Array{Int,1}(), Array{Int,2}() ))
    Q,R,E = qr(X, Val(true))
    diagr = abs.(diag(R))
    r = findall(diagr .>= tol*diagr[1])[end]
    idx = sort(E[1:r])
    return idx
end


function echelonize(X::T,v; eps::Real=1e-10) where T <:DenseArray
    M,N = size(X)

    idxrow = idxlicols(Matrix(X'))
    Mred = length(idxrow)

    idxind = idxlicols(X)
    idxdep = setdiff(1:N,idxind)
    newidx = vcat(idxind,idxdep)
    Tv = @view X[idxrow,newidx] 
    iTv = inv(Tv[1:Mred, 1:Mred])
    res = iTv * Tv
    # trimming zeros
    for i in eachindex(res)
        abs(res[i]) < eps && (res[i] = zero(res[i]))
    end
    bnew = iTv * v[idxrow]
    # trimming ones
    for i in 1:Mred
        abs(1.0 - res[i,i]) < eps && (res[i,i] = one(res[i,i]))
    end

    return idxdep, newidx, res, bnew

end

function preprocess(S,b,lb,ub,rxns)
    b,lb, ub = copy(b), copy(lb), copy(ub)
    n = length(lb)
    ei = zeros(n)
    for i=1:n
        ei[i] = -1.0
        sol = linprog(ei, S, b, b, lb, ub, ClpSolver())
        ub[i] = min(ub[i], sol.sol[i])
        ei[i] = +1.0
        sol = linprog(ei, S, b, b, lb, ub, ClpSolver())
        lb[i] = max(lb[i], sol.sol[i])
        ei[i] = 0.0
    end
    blocked = (lb .== ub)
    idxb = findall(lb .== ub)
    println("$(sum(blocked)) blocked fluxes")
    for i = 1:sum(blocked)
        println(rxns[idxb[i]], "is fixed to ", lb[idxb[i]])
    end
    unblocked = (lb .< ub)
    return S[:,unblocked], b - S*(blocked .* lb), lb[unblocked], ub[unblocked], rxns[unblocked]
end
