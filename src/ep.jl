Φ(x) = 0.5*(1+erf(x/sqrt(2.)))
ϕ(x) = exp(-x.^2/2)/sqrt(2π)

function inplaceinverse!(dest::AbstractArray,source::AbstractArray)
    dest = copy!(dest, source)
    Base.LinAlg.inv!(cholfact!(Hermitian(dest)))
end

# K * x == Y 
# K stoichiometrix matrix M x N (M metabolites, N fluxes)
# Y = N dimensional vector
# nusup, nuinf =  N-dimensional vectors with fluxes' lower and upper bounds
metabolicEP(X::COBRA.LPproblem; args...) = metabolicEP(X.S, X.b, X.lb, X.ub; args...) # convinence for runnning directly on Lpproblems

function metabolicEP{T<:AbstractFloat}(K::AbstractArray{T,2}, Y::Array{T,1}, nuinf::Array{T,1}, nusup::Array{T,1};
                                       beta::Real=1e7,      # inverse temperature
                                       verbose::Bool=true,  # output verbosity
                                       damp::Real=0.9,      # damp ∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield
                                       epsconv::Real=1e-6,  # convergence criterion
                                       maxiter::Int=2000,   # maximum iteration count
                                       maxvar::Real=1e50,   # maximum numerical variance
                                       minvar::Real=1e-50,  # minimum numerical variance
                                       solution::Union{EPout{T},Void}=nothing,  # start from a solution
                                       expval=nothing)      # fix posterior probability experimental values for std and mean
   

    updatefunction,scalefact,epfield = prepareinput(K,Y,nuinf,nusup,beta,verbose,solution,expval,T)
    
    epalg = EPAlg(beta, minvar, maxvar, epsconv, damp, maxiter,verbose)    
    updatealg = beta < Inf ? eponesweep! : eponesweepT0!
    epmat = beta < Inf ? EPMat(K,Y,nuinf, nusup, beta) : EPMatT0(K,Y,nuinf, nusup)
    
    
    returnstatus=epconverge!(epfield,epmat,epalg, updatealg)
    scaleepfield!(epfield,nusup,nuinf,Y,scalefact)
    return  EPout(epfield.μ,epfield.s, epfield.av, epfield.va, epfield, returnstatus)
end


function prepareinput(K,Y,nuinf,nusup,beta,verbose,solution,expval,T)

    M,N = size(K) 
    M < N || warn("M = $M ≥ N = $N")
    sum(nusup .< nuinf) == 0 || error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")     

    verbose && println("Analyzing a $M x $N stoichiometric matrix.")
    
    updatefunction = if beta == Inf
        if isstandardform(K)
            eponesweepT0!
        else
            error("for T = 0 algorithm, S should be [I | M] whre I is the Identity and M any matrix")
        end
    else
        eponesweep!
    end

    scalefact = max(maximum(abs.(nuinf)), maximum(abs.(nusup)))
    if solution === nothing
        epfield = EPFields(N,expval,scalefact,T)
    else
        epfield = solution.sol
    end


    scale!(nusup,1.0/scalefact)
    scale!(nuinf,1.0/scalefact)
    scale!(Y,1.0/scalefact)
    
    return updatefunction,scalefact,epfield
end

function epconverge!{T<:Function,M<:AbstractEPMat}(epfield::EPFields,epmat::M,epalg::EPAlg, eponesweep!::T)

    @extract epalg : maxiter verbose epsconv
    
    returnstatus = :unconverged
    iter = 0
    while iter < maxiter
        iter += 1
        (errav,errvar,errμ, errs) = eponesweep!(epfield,epalg, epmat)
        verbose && @printf("it = %d beta = %g errav = %g errvar = %g errμ = %g errs = %g\n",
                           iter, epalg.beta, errav, errvar, errμ, errs)
        if max(errav, errvar,errμ,errs) < epsconv
            returnstatus = :converged
            break
        end
    end
    return returnstatus
end



function scaleepfield!(X,nuinf,nusup,Y,scalefact)
    @extract X : μ s av va
    scale!(μ, scalefact)
    scale!(s, scalefact^2)
    scale!(av, scalefact)
    scale!(va, scalefact^2)
    scale!(nusup,scalefact)
    scale!(nuinf,scalefact)
    scale!(Y,scalefact)
end

function eponesweepT0!(X::EPFields, epalg::EPAlg, epmatT0::EPMatT0)
    @extract X : av va a b μ s siteflagave siteflagvar
    @extract epalg : beta minvar maxvar epsconv damp
    @extract epmatT0 : Σy Σw G nuinf nusup  
    
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
    errav,errva,errμ,errs = (minerr,minerr,minerr,minerr)

    Σw = inv(Diagonal(1.0 ./ bw) + G' * Diagonal( 1.0 ./ by ) * G)
    for i in 1:N-M vaw[i] = Σw[i,i] end
    
    A_mul_B!(Σy,G*Σw,G')
    for i in 1:M vay[i] = Σy[i,i] end
    
    A_mul_B!(avw,Σw, aw ./ bw + G'*(ay ./ by))
    A_mul_B!(avy,G,μw)

    for i in eachindex(μw)  # loop M+1:N
        newμw,newsw = newμs(vaw[i],aw[i],bw[i],avw[i],nuinf[i+M],nusup[i+M], minvar,maxvar)
        errμ = max(errμ, abs(μw[i]-newμw))
        errs = max(errs, abs(sw[i]-newsw))
        μw[i] = newμw
        sw[i] = newsw
        
        newavw,newvaw = newav(sw[i],μw[i],avw[i],vaw[i],siteflagave[i+M],siteflagvar[i+M],
                              nuinf[i+M],nusup[i+M],minvar,maxvar)
        errav = max(errav,abs(avw[i]-newavw))
        errva = max(errva,abs(vaw[i]-newvaw))
        avw[i] = newavw
        vaw[i] = newvaw 

        newaw,newbw = matchmom(μw[i],sw[i],avw[i],vaw[i],minvar,maxvar)
        aw[i] = damp * aw[i] + (1.0-damp)*newaw
        bw[i] = damp * bw[i] + (1.0-damp)*newbw
    end
    
    for i in eachindex(μy)   # loop  1:M
        newμy,newsy = newμs(vay[i],ay[i],by[i],avy[i],nuinf[i],nusup[i],minvar,maxvar)
        errμ = max(errμ, abs(μy[i]-newμy))
        errs = max(errs, abs(sy[i]-newsy))
        μy[i] = newμy
        sy[i] = newsy

        newavy,newvay = newav(sy[i],μy[i],avy[i],vay[i],siteflagave[i],siteflagvar[i],
                              nuinf[i],nusup[i],minvar,maxvar)
        errav = max(errav,abs(avy[i]-newavy))
        errva = max(errva,abs(vay[i]-newvay))
        avy[i] = newavy
        vay[i] = newvay 
        
        neway,newby = matchmom(μy[i],sy[i],avy[i],vay[i],minvar,maxvar)
        ay[i] = damp * ay[i] + (1.0-damp)*neway
        by[i] = damp * by[i] + (1.0-damp)*newby
    end    
    return errav,errav, errμ, errs
end    

function matchmom(μ,s,av,va, minvar,maxvar)
    newb = clamp(inv(1.0/va - 1.0/s),minvar,maxvar)
    newa = av + newb*(av-μ)/s
    isnan(newa) || isnan(newb) && warn("a = $newa b = $newb")
    return newa, newb
end

function newav(s,μ,av,va,siteflagave,siteflagvar,nuinf,nusup, minvar, maxvar)
    sqrts = sqrt(s)
    xinf = (nuinf - μ) / sqrts
    xsup = (nusup - μ) / sqrts
    scra1,scra12 = compute_mom5d(xinf,xsup)
    avnew  = siteflagave ? μ + scra1*sqrts : av 
    varnew = siteflagvar ? max(minvar,s*(1.0+scra12)) : va    
    isnan(avnew) || isnan(varnew) && println("avnew = $avnew varnew = $varnew")
    return avnew, varnew
end

function newμs(Σ,a,b,v,nuinf,nusup,minvar,maxvar)
    s = clamp(inv(1.0/Σ - 1.0/b),minvar,maxvar)   
    μ = if Σ != b
        s * (v/Σ - a/b)
    else
        warn("I'm here: nusup = ",nusup," nuinf = ",nuinf, " Σ = ", Σ)
        0.5 * (nusup+nuinf)
    end
    return μ,s
end

let DDwXDy = Dict{Int,Matrix}()
    global fast_similarity_inv!
    function fast_similarity_inv!{T<:AbstractFloat}(dest,Dw::Vector{T},Dy::Vector{T},G)        
        NmM = length(Dw)
        DwXDy = Base.get!(DDwXDy,NmM,zeros(T,NmM,NmM))
        @inbounds for i in eachindex(Dw)
            DwXDy[i,i] = 1.0 / Dw[i]
        end
        BLAS.syrk!('U','T',1.0, Diagonal((1./sqrt.(Dy)))*G,1.0,DwXDy)
        
        inplaceinverse!(dest, DwXDy)
        return nothing
    end
end

function eponesweep!(X::EPFields,epalg::EPAlg, epmat::EPMat)
    @extract X : av va a b μ s siteflagave siteflagvar
    @extract epalg : beta minvar maxvar epsconv damp
    @extract epmat : KK KKPD invKKPD nuinf nusup D KY v
    @extract epmat : invKKPD KKPD KK D
    
    for i in eachindex(D)
        KKPD[i,i] = KK[i,i] + D[i]
    end
    inplaceinverse!(invKKPD,KKPD)
    
    T = eltype(v)
    
    errav = typemin(T)
    errva = typemin(T)
    errμ  = typemin(T)
    errs  = typemin(T)

    A_mul_B!(v,invKKPD, (KY + D.*a))
    
    for i in eachindex(av)
        newμ,news = newμs(invKKPD[i,i],a[i],b[i],v[i],nuinf[i],nusup[i],minvar, maxvar)
        errμ = max(errμ, abs(μ[i]-newμ))
        errs = max(errs, abs(s[i]-news))
        μ[i] = newμ
        s[i] = news
        
        newave,newva = newav(s[i],μ[i],av[i],va[i],siteflagave[i],siteflagvar[i],nuinf[i],nusup[i],minvar,maxvar)
        errav = max(errav,abs(av[i]-newave))
        errva = max(errva,abs(va[i]-newva))
        av[i] = newave
        va[i] = newva
        
        newa,newb = matchmom(μ[i],s[i],av[i],va[i],minvar,maxvar)
        a[i] = damp * a[i] + (1.0-damp)*newa
        b[i] = damp * b[i] + (1.0-damp)*newb
        D[i] = 1.0/b[i]
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

function _parseexpval!{T<:Tuple}(expval::T,siteflagave::BitArray{1},siteflagvar::BitArray{1})

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

function _parseexpval!{T<:Vector}(expval::T,siteflagave::BitArray{1},siteflagvar::BitArray{1})

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



