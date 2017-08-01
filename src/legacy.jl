
function eponesweepold!(X::EPFields, epalg::EPAlg, epmat::EPMat)
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

