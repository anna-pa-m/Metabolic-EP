immutable EPFields{T<:AbstractFloat}
    av::Vector{T}
    va::Vector{T}
    a::Vector{T}
    b::Vector{T}
    μ::Vector{T}
    s::Vector{T}
    siteflagave::BitArray{1}
    siteflagvar::BitArray{1}
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
             zeros(T,N),
             ones(T,N),
             siteflagave,
             siteflagvar)
end


function EPFields(N::Int,expval::Void,scalefact,T)    

    return EPFields(zeros(T,N),
             zeros(T,N),
             zeros(T,N),
             ones(T,N),
             zeros(T,N),
             ones(T,N),
             trues(N),
             trues(N))                
end

abstract type AbstractEPMat end

struct EPMat{T<:AbstractFloat} <: AbstractEPMat
    KK::AbstractArray{T,2}
    KKPD::AbstractArray{T,2}
    invKKPD::Matrix{T}
    KY::Vector{T}
    v::Vector{T}
    lb::Vector{T}
    ub::Vector{T}
end

function EPMat{T<:AbstractFloat}(K::AbstractArray{T}, Y::Vector{T}, lb::Vector{T}, ub::Vector{T}, beta::T)
    M,N = size(K)
    KKPD = full(beta * K' * K)
    if beta != Inf
        return EPMat(copy(KKPD), copy(KKPD), zeros(T,N,N), beta * K' * Y, zeros(T,N),lb,ub)
    else
        error("I really should not be here")
    end
end


struct EPMatT0{T<:AbstractFloat} <: AbstractEPMat
    Σy::Matrix{T}
    Σw::Matrix{T}
    G::Matrix{T}
    lb::Vector{T}
    ub::Vector{T}
    vy::Vector{T}
    vw::Vector{T}
    Y::Vector{T}
end

function EPMatT0{T<:AbstractFloat}(K::AbstractArray{T,2}, Y::Vector{T}, lb::Vector{T}, ub::Vector{T})
    M,N = size(K)
    M <= N || error("M=$M cannot be larger than N=$N")
    return EPMatT0(zeros(T,M,M), zeros(T,N-M,N-M), copy(full(K[1:M,M+1:N])), lb,ub,zeros(T,M),zeros(T,N-M),Y)
end

mutable struct EPAlg{T<:AbstractFloat}
    beta::T
    minvar::T
    maxvar::T
    epsconv::T
    damp::T
    maxiter::Int
    verbose::Bool
end

struct EPout{T<:AbstractFloat}
    μ::Vector{T}
    σ::Vector{T}
    av::Vector{T}
    va::Vector{T}
    sol::EPFields{T}
    status::Symbol
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
