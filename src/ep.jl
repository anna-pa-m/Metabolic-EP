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
   
    M,N = size(K) 
    M < N || warn("M = $M ≥ N = $N")
    sum(nusup .< nuinf) == 0 || error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")     
    verbose && println("Analyzing a $M x $N stoichiometric matrix.")

    
    updatefunction = if beta == Inf
        if isstandardform(K)
            eponesweepT0!
        else
            error("for T = 0 algorithm, S should be [I | M] whre I is diagonal")
        end
    else
        eponesweep!
    end

    
    siteflagave = trues(N)
    siteflagvar = trues(N)    
    scalefact = max(maximum(abs.(nuinf)), maximum(abs.(nusup)))

    scale!(nusup,1.0/scalefact)
    scale!(nuinf,1.0/scalefact)
    scale!(Y,1.0/scalefact)
    
    
    if solution === nothing
        epfield = EPFields(N,expval,scalefact,T)
    else
        epfield = solution.sol
    end

    epalg = EPAlg(beta, minvar, maxvar, epsconv, damp, maxiter,verbose)    
    epmat = EPMat(K,Y,nuinf, nusup, beta)

    returnstatus=epconverge!(epfield,epmat,epalg, eponesweep!)    

    scaleepfield!(epfield,scalefact)
    scale!(nusup,scalefact)
    scale!(nuinf,scalefact)
    scale!(Y,scalefact)

    return  EPout(epfield.μ,epfield.s, epfield.av, epfield.va, epfield, returnstatus)
end

function epconverge!{T<:Function}(epfield::EPFields,epmat::EPMat,epalg::EPAlg, eponesweep!::T)

    @extract epalg : maxiter verbose epsconv
    
    returnstatus = :unconverged
    iter = 0
    while iter < maxiter
        iter += 1
        updatemat!(epmat)        
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

function scaleepfield!(X::EPFields,scalefact)
    @extract X : μ s av va
    scale!(μ, scalefact)
    scale!(s, scalefact^2)
    scale!(av, scalefact)
    scale!(va, scalefact^2)
end

function updatemat!(epmat)
    @extract epmat : invKKPD KKPD KK D I
        
    for i in eachindex(D)
        KKPD[i,i] = KK[i,i] + D[i]
    end
    
    inplaceinverse!(invKKPD,KKPD)

    for i in eachindex(D)
        I[i] = invKKPD[i,i]
    end
   
    nothing
end



function eponesweepT0!(X::EPFields, epalg::EPAlg, epmat::EPMat)
    @extract X : av va a b μ s  siteflagave siteflagvar
    @extract epalg : beta minvar maxvar epsconv damp
    @extract epmat : KK KKPD invKKPD nuinf nusup D KY I v

    





end


function eponesweep!(X::EPFields, epalg::EPAlg, epmat::EPMat)
    @extract X : av va a b μ s  siteflagave siteflagvar
    @extract epalg : beta minvar maxvar epsconv damp
    @extract epmat : KK KKPD invKKPD nuinf nusup D KY I v

    T = eltype(v)
    
    errav = typemin(T)
    errvar = typemin(T)
    errμ  = typemin(T)
    errs  = typemin(T)

    A_mul_B!(v,invKKPD, (KY + D.*a))
    
    for i in eachindex(I)
        I[i] = min(I[i],b[i])
        s1   = max(minvar,1/I[i] - 1/b[i])
        news = max(minvar,1/s1) 
        errs = max(errs,abs(news-s[i])) 
        s[i] = news
        if I[i] != b[i]
            newμ = (v[i]-a[i]*I[i]/b[i])/(1.0-I[i]/b[i])
            errμ = max(errμ,abs(μ[i]-newμ))
            isinf(errμ) && (println("μ[$i] = ", μ[i]," newμ = $newμ"); break) 
            μ[i] = newμ
        else
            warn("I'm here: nusup[$i] = ",nusup[i]," nuinf[$i] = ",nuinf[i], " I[$i] = ",I[i])
            newμ = 0.5 * (nusup[i] + nuinf[i])
            errμ = max(errμ,abs(μ[i]-newμ))
            μ[i] = newμ
        end

        sqrts  = sqrt(s[i])
        xsup   = (nusup[i] - μ[i])/sqrts
        xinf   = (nuinf[i] - μ[i])/sqrts

        scra1,scra12 = compute_mom5d(xinf,xsup)
        if siteflagave[i]
            avnew = μ[i] + scra1 * sqrts
            errav = max(errav, abs(av[i] - avnew))               
        else
            avnew = av[i]
        end
        if siteflagvar[i]
            varnew = max(minvar, s[i] * (1.0 + scra12))                                        
            errvar = max(errvar, abs(va[i]-varnew))
        else
            varnew = va[i]
        end

        av[i]  = avnew
        va[i] = varnew            
        
        isnan(avnew) || isnan(varnew) && println("avnew = $avnew varnew = $varnew")
        isnan(a[i]) || isnan(b[i])  && println("a[$i] = ", a[i]," b[$i] = ",b[i])
            
        new_b = min(maxvar, max(minvar, 1./(1./va[i] - 1./s[i])))
        new_a = av[i] + new_b*(av[i] - μ[i])/s[i]
        a[i] = damp*a[i]  + (1.0 - damp)*new_a
        b[i] = damp*b[i]  + (1.0 - damp)*new_b            
    end
    
    for i in eachindex(D)
        D[i] = 1./b[i]
    end
    
    return (errav,errvar,errμ,errs)
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



