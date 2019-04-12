using LinearAlgebra
Φ(x) = 0.5*(1.0+erf(x/sqrt(2.0)))
ϕ(x) = exp(-x.^2/2.0)/sqrt(2π)

function inplaceinverse!(dest::AbstractArray,source::AbstractArray)
    dest = copyto!(dest, source)
    LinearAlgebra.inv!(cholesky!(Hermitian(dest)))
end


"""
res=metabolicEP(S,b,lb,ub,...)


The output in res is of type `EPout`: there are several fields:
-   `μ::Vector`: A parameter linked to the mean of the posterior probability 
-   `σ::Vector`: A parameter linked to the std  of the posterior probability 
-   `av::Vector`: The mean posterior probability
-   `va::Vector`: The variance of the posterior probability
-   `sol::EPFields`: The internal field status. From this value we can restart the sampling from a specific state.
-   `status::Symbol`: either ``:converged`` or ``:unconverged``.

Input (required)
----
- `S`: MxN matrix (either sparse or dense) please note that if you input a dense version, the algorithm is slighlty more efficient. Dense matrices can be create from sparse ones with `Matrix(S)`.
- `b`: a vector of M intakes/uptakes 
- `lb`: a vector of lengh N of lower bounds.
- `ub`: a vector of lengh N of upper bounds.

Input (optional arguments). 
----
- `beta` (inverse temperature::``Real``): default 10^7 
- `verbose` (``true`` or ``false``): default ``true``
- `damp` (∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield): default 0.9  
- `epsconv` (convergence criterion): default 1e-6
- `maxiter` (maximum number of iterations): default 2000
- `maxvar`  (threshold on maximum variance): default 1e50
- `minvar`  (threshold on minimum variance): default 1e-50
- `solution` (start from solution. Is of type ``EPout``): default: ``nothing``
- `expval` (fix to posterior probability of mean and/or variance to values): default ``nothing``. expval can be either at ``Tuple{Float64,Float64,Int}`` or a ``Vector{Tuple{Float64,Float64,Int}}``. Values can be fixed as``expval=(0.2,0.4,4)`` meaning that for flux index 4 the mean is set to 0.2 and the variance to 0.4. Fixing more values ``expval=[(0.2, 0.3, 4), (0.4, nothing, 5)]``: in this case, we fix the posterior of flux 4 to 0.2 (mean) and 0.3 (variance), while for flux 5 we fix the mean to 0.4 and we keep the variance free.

"""
function metabolicEP(K::AbstractArray{T,2}, Y::Array{T,1}, lb::Array{T,1}, ub::Array{T,1};
                     beta::Real=1e7,      # inverse temperature
                     verbose::Bool=true,  # output verbosity
                     damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
                     epsconv::Real=1e-6,  # convergence criterion
                     maxiter::Int=2000,   # maximum iteration count
                     maxvar::Real=1e50,   # maximum numerical variance
                     minvar::Real=1e-50,  # minimum numerical variance
                     solution::Union{EPout{T},Nothing}=nothing,  # start from a solution
                     expval=nothing # fix posterior probability experimental values for std and mean
                     ) where T<:Real


    llb = copy(lb) # making  a local copy to rescale
    lub = copy(ub)
    updatealg,scalefact,epfield = prepareinput(K,Y,llb,lub,beta,verbose,solution,expval,T)
    scaleepfield!(epfield,lub,llb,Y,1.0/scalefact) # rescaling fields in [0,1]
    epalg = EPAlg(beta, minvar, maxvar, epsconv, damp, maxiter,verbose)
    epmat = beta < Inf ? EPMat(K,Y,llb, lub, beta) : EPMatT0(K,Y,llb, lub)
    returnstatus = epconverge!(epfield,epmat,epalg, updatealg)
    scaleepfield!(epfield,lub,llb,Y,scalefact,beta)
    if beta < Inf
        return  EPout(epfield.μ,epfield.s, epfield.av, epfield.va, epfield, returnstatus)
    else
        idx = epmat.idx
        return  EPout(epfield.μ[idx],epfield.s[idx], epfield.av[idx], epfield.va[idx], epfield, returnstatus)
    end
end

function prepareinput(K,Y,lb,ub,beta,verbose,solution,expval,T)

    M,N = size(K)
    M < N || @warn("M = $M ≥ N = $N")
    all(ub .< lb) || error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds") 

    verbose && println(stderr, "Analyzing a $M x $N stoichiometric matrix.")

    scalefact = T(0)
    
    updatefunction = if beta == Inf
        eponesweepT0!
    else
        eponesweep!
    end

    scalefact = max(maximum(abs.(lb)), maximum(abs.(ub)))
    if solution === nothing
        epfield = EPFields(N,expval,scalefact,T)
    else
        epfield = deepcopy(solution.sol) # preserve the original solution!
    end

    return updatefunction,scalefact,epfield
end

function epconverge!(epfield::EPFields,epmat::M,epalg::EPAlg, eponesweep!::T) where {T<:Function,M<:AbstractEPMat}


    @extract epalg : maxiter verbose epsconv
    
    returnstatus = :unconverged
    iter = 0
    print(stderr, "Converging with β=$(epalg.beta) maxth=$(epsconv) maxiter=$(maxiter):\n")
    while iter < maxiter
        iter += 1
        (errav,errvar,errμ, errs) = eponesweep!(epfield,epalg, epmat)
        if (verbose)
            @printf(stderr, "it = %d ɛav = %.2e ɛvar = %.2e ɛμ = %.2e ɛs = %.2e                 \r", iter, errav, errvar, errμ, errs)
            flush(stderr)
        end
        
        if max(errav, errvar) < epsconv
            returnstatus = :converged
            break
        end
    end
    if verbose
        print(stderr, "\n")
        flush(stderr)
    end    
    return returnstatus
end



function scaleepfield!(epfield,ub,Y,scalefact,beta)

    if beta < Inf
        @extract epfield : μ s av va
        idx = [i for i ∈ 1:length(av) ]        
    else
        @extract epfield : μ s av va idx
    end

    rmul!(μ, scalefact)
    rmul!(s, scalefact^2)
    rmul!(av, scalefact)
    rmul!(va, scalefact^2)
    rmul!(ub,scalefact)
    rmul!(lb,scalefact)
    rmul!(Y,scalefact)
end

function eponesweepT0!(epfields::EPFields, epalg::EPAlg, epmatT0::EPMatT0)
    @extract epfields : av va a b μ s siteflagave siteflagvar
    @extract epalg : beta minvar maxvar epsconv damp
    @extract epmatT0 : Σy Σw G lb ub vy vw Y 

    M = size(G,1)
    N = length(av)
    
    idxy = 1:M
    idxw = M+1:N

    ay,aw = view(a,idxy), view(a,idxw)
    by,bw = view(b,idxy), view(b,idxw)
    sy,sw = view(s,idxy), view(s,idxw)
    μy,μw = view(μ,idxy), view(μ,idxw)
    avy,avw = view(av,idxy), view(av,idxw)
    vay,vaw = view(va,idxy), view(va,idxw)

    minerr = typemin(av[1])
    errav,errva,errμ,errs = minerr,minerr,minerr,minerr

    Σw = inv(Diagonal(1.0 ./ bw) + G' * Diagonal( 1.0 ./ by ) * G)
    #fast_similarity_inv!(Σw, bw,  by, G)
    mul!(Σy,G*Σw,G')
    mul!(vw,Σw, aw ./ bw - G'*(ay ./ by))
    mul!(vy,G,vw)
    for i in eachindex(vy) vy[i] = -vy[i] + Y[i] end

    for i in eachindex(μw)  # loop M+1:N
        newμw,newsw = newμs(Σw[i,i],aw[i],bw[i],vw[i],lb[i+M],ub[i+M], minvar,maxvar)
        errμ = max(errμ, abs(μw[i]-newμw))
        errs = max(errs, abs(sw[i]-newsw))
        μw[i] = newμw
        sw[i] = newsw
        # println("μw[$(i+M)] = ", μw[i]," sw[$(i+M)] = ", sw[i], " Σw = ",Σw[i,i] )

        
        newavw,newvaw = newav(sw[i],μw[i],avw[i],vaw[i],siteflagave[i+M],siteflagvar[i+M],
                              lb[i+M],ub[i+M],minvar,maxvar)
        errav = max(errav,abs(avw[i]-newavw))
        errva = max(errva,abs(vaw[i]-newvaw))
        avw[i] = newavw
        vaw[i] = newvaw 

        newaw,newbw = matchmom(μw[i],sw[i],avw[i],vaw[i],minvar,maxvar)
        aw[i] = damp * aw[i] + (1.0-damp)*newaw
        bw[i] = damp * bw[i] + (1.0-damp)*newbw
    end
    
    for i in eachindex(μy)   # loop  1:M
        
        newμy,newsy = newμs(Σy[i,i],ay[i],by[i], vy[i],lb[i],ub[i],minvar,maxvar)
        errμ = max(errμ, abs(μy[i]-newμy))
        errs = max(errs, abs(sy[i]-newsy))
        μy[i] = newμy
        sy[i] = newsy
#        println("μy[$i] = ", μy[i]," sy[$i] = ", sy[i], " Σ = ", Σy[i,i], " (",lb[i],":",ub[i],")"," ay[$i] = ",ay[i], " by[$i] = ", by[i])
        
        newavy,newvay = newav(sy[i],μy[i],avy[i],vay[i],siteflagave[i],siteflagvar[i],
                              lb[i],ub[i],minvar,maxvar)
        errav = max(errav,abs(avy[i]-newavy))
        errva = max(errva,abs(vay[i]-newvay))
        avy[i] = newavy
        vay[i] = newvay 
        
        neway,newby = matchmom(μy[i],sy[i],avy[i],vay[i],minvar,maxvar)
        ay[i] = damp * ay[i] + (1.0-damp)*neway
        by[i] = damp * by[i] + (1.0-damp)*newby
    end    
    return errav, errva, errμ, errs
end

function matchmom(μ,s,av,va, minvar,maxvar)
    newb = clamp(inv(1.0/va - 1.0/s),minvar,maxvar)
    newa = av + newb*(av-μ)/s
    isnan(newa) || isnan(newb) && @warn("a = $newa b = $newb")
    return newa, newb
end

function newav(s,μ,av,va,siteflagave,siteflagvar,lb,ub, minvar, maxvar)
    sqrts = sqrt(s)
    xinf = (lb - μ) / sqrts
    xsup = (ub - μ) / sqrts
    scra1,scra12 = compute_mom5d(xinf,xsup)
    avnew  = siteflagave ? μ + scra1*sqrts : av 
    varnew = siteflagvar ? max(minvar,s*(1.0+scra12)) : va    
    isnan(avnew) || isnan(varnew) && println("avnew = $avnew varnew = $varnew")
    return avnew, varnew
end

function newμs(Σ,a,b,v,lb,ub,minvar,maxvar)

    Σ == 0 && (Σ = minvar)
    #lΣ = clamp(Σ,minvar,maxvar)
    #s = Σ > 0 ? clamp(inv(1.0/Σ - 1.0/b),minvar,maxvar) : minvar   
    s = clamp(inv(1.0/Σ - 1.0/b),minvar,maxvar)   
    μ = if Σ != b 
        s * (v/Σ - a/b)
    else
        #@warn("I'm here: ub = ",ub," lb = ",lb, " Σ = ", Σ)
        0.5 * (ub+lb)
    end
    return μ,s
end

let DDwXDy = Dict{Int,Matrix}()
    global fast_similarity_inv!
    function fast_similarity_inv!(dest::Matrix{T},Dw,Dy,G) where T <: Real
        NmM = length(Dw)
        DwXDy = Base.get!(DDwXDy,NmM,zeros(T,NmM,NmM))
        fill!(DwXDy,zero(T))
        @inbounds for i in eachindex(Dw)
            DwXDy[i,i] = 1.0 / Dw[i]
        end
        BLAS.syrk!('U','T',1.0, Diagonal((1.0/sqrt.(Dy)))*G,1.0,DwXDy)        
        inplaceinverse!(dest, DwXDy)
        return nothing
    end
end


function eponesweep!(X::EPFields,epalg::EPAlg, epmat::EPMat)
    @extract X : av va a b μ s siteflagave siteflagvar
    @extract epalg : beta minvar maxvar epsconv damp
    @extract epmat : KK KKPD invKKPD lb ub KY v
    
    for i in eachindex(b)
        KKPD[i,i] = KK[i,i] + 1.0/b[i]
    end
    inplaceinverse!(invKKPD,KKPD)

    minerr = typemin(av[1])
    errav,errva,errμ,errs = minerr,minerr,minerr,minerr

    mul!(v,invKKPD, (KY + a./b))
    
    for i in eachindex(av)
        newμ,news = newμs(invKKPD[i,i],a[i],b[i],v[i],lb[i],ub[i],minvar, maxvar)
        errμ = max(errμ, abs(μ[i]-newμ))
        errs = max(errs, abs(s[i]-news))
        μ[i] = newμ
        s[i] = news


        newave,newva = newav(s[i],μ[i],av[i],va[i],siteflagave[i],siteflagvar[i],lb[i],ub[i],minvar,maxvar)
        errav = max(errav,abs(av[i]-newave))
        errva = max(errva,abs(va[i]-newva))
        av[i] = newave
        va[i] = newva
        
        newa,newb = matchmom(μ[i],s[i],av[i],va[i],minvar,maxvar)
        a[i] = damp * a[i] + (1.0-damp)*newa
        b[i] = damp * b[i] + (1.0-damp)*newb
    end
    return errav,errva,errμ,errs
end


function compute_mom5d(xinf, xsup)
  
    minval = min(abs(xinf), abs(xsup))
    sgn = sign(xinf*xsup)

    if xsup - xinf < 1e-8
        return (0.5*(xsup + xinf),-1.)
    end

    if minval <= 6. || sgn <= 0
        ϕsup   = ϕ(xsup)
        Φsup   = Φ(xsup)
        ϕinf   = ϕ(xinf)
        Φinf   = Φ(xinf)
        scra1 = (ϕinf - ϕsup)/(Φsup - Φinf)             
        scra2 = (xinf * ϕinf - xsup*ϕsup)/(Φsup -Φinf)
        scra12 = scra2 - scra1^2
        return scra1, scra12
    else
        delta2 = (xsup^2 - xinf^2)*0.5 
        if delta2 > 40.
            scra1 = xinf^5/(3 - xinf^2 + xinf^4)
            scra2 = xinf^6/(3 - xinf^2 + xinf^4)
        else
            scra1 = (xinf*xsup)^5 * (1. - exp(delta2)) / (-exp(delta2)*(3.0-xinf^2 + xinf^4)*xsup^5 + xinf^5*(3-xsup^2 + xsup^4))
            scra2 = (xinf*xsup)^5 * (xsup - xinf*exp(delta2)) / (-exp(delta2)*(3.0-xinf^2 + xinf^4)*xsup^5 + xinf^5*(3-xsup^2 + xsup^4))
        end
        scra12 = scra2 - scra1^2        
        isnan(scra1) || isnan(scra2) || isnan(scra12) && println("scra1 = $scra1 scra2 = $scra2")
        !isfinite(scra1) ||  !isfinite(scra2) && println("scra1 = $scra1 scra2 = $scra2")
        return scra1, scra12
    end
end


function parseexpval!(expval,siteflagave::BitArray{1}, siteflagvar::BitArray{1},scalefact::Float64)

    expave,expvar=_parseexpval!(expval,siteflagave,siteflagvar)

    for (k,v) in expave
        expave[k] = v/scalefact
    end

    for (k,v) in expvar
        expvar[k] = v/scalefact^2
    end
    expave,expvar
end

function _parseexpval!(expval::Tuple,siteflagave::BitArray{1},siteflagvar::BitArray{1})

    N = length(siteflagave)
    length(expval) == 3 || error("We expect a 3-uple here")
    expsite = expval[3]
    1<= expsite <= N || error("expsite = $expsite not ∈ 1,...,$N")
    expave = Dict{Int,Float64}()
    expvar = Dict{Int,Float64}()    

    if expval[1] != nothing
        siteflagave[expsite] = false
        expave[expsite] = expval[1]
    end
    if expval[2] != nothing
        siteflagvar[expsite] = false
        expvar[expsite] = expval[2]
    end
    expave,expvar
end

function _parseexpval!(expval::Vector,siteflagave::BitArray{1},siteflagvar::BitArray{1})

    N = length(siteflagave)    
    expave = Dict{Int,Float64}()
    expvar = Dict{Int,Float64}()
    for i in eachindex(expval)
        expsite = expval[i][3]
        1 <= expsite <= N || error("expsite = $minsite not ∈ 1,...,$N")
        if expval[i][1] != nothing
            siteflagave[expsite] = false
            expave[expsite] = expval[i][1]
        end
        if expval[i][2] != nothing
            siteflagvar[expsite] = false
            expvar[expsite] =  expval[i][2]
        end
    end

    expave,expvar
end

_parseexpval!(nothing,siteflagave::BitArray{1},siteflagvar::BitArray{1})=(Dict{Int,Float64}(),Dict{Int,Float64}())
