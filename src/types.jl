immutable EPFields{T<:AbstractFloat}
    av::Vector{T}
    var::Vector{T}
    a::Vector{T}
    b::Vector{T}
    D::Vector{T}
    μ::Vector{T}
    s::Vector{T}
    new_a::Vector{T}
    new_b::Vector{T}
    siteflagave::BitArray{1}
    siteflagvar::BitArray{1}
end


type EPout{T<:AbstractFloat}
    μ::Vector{T}
    σ::Vector{T}
    av::Vector{T}
    va::Vector{T}
    sol::EPFields{T}
    status::Symbol
end


function EPFields(N::Int,expval,scalefact,T)

    siteflagvar = trues(N)
    siteflagave = trues(N)
    
    expave, expvar = parseexpval!(expval,siteflagave,siteflagvar,scalefact)
    av = zeros(T,N)
    var = zeros(T,N)

    for (k,v) in expave
        av[k] = v
    end    
    for (k,v) in expvar
        var[k] = v
    end

    EPFields(av,
             var,
             zeros(T,N),
             ones(T,N),
             ones(T,N),
             zeros(T,N),
             ones(T,N),
             zeros(T,N),
             zeros(T,N),
             siteflagave,
             siteflagvar)
end


function EPFields(N::Int,expval::Void,scalefact,T)    

    EPFields(zeros(T,N),
             zeros(T,N),
             zeros(T,N),
             ones(T,N),
             ones(T,N),
             zeros(T,N),
             ones(T,N),
             zeros(T,N),
             zeros(T,N),
             trues(N),
             trues(N))            
    
end


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
