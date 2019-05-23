module TestEP
using MetabolicEP, Test

function testEP()
    S = [1.0 1.0 1.0 -1.0]
    b=zeros(size(S,1))
    lb=zeros(size(S,2))
    ub=ones(size(S,2));
    truemean = [1//4, 1//4, 1//4, 3//4]
    res0=metabolicEP(S,b,lb,ub,beta=Inf,epsconv=1e-9,verbose=false)
    res=metabolicEP(S,b,lb,ub,beta=1e9,epsconv=1e-9,verbose=false)
    @test res0.status == :converged
    @test res.status == :converged
    @test res.av ≈ truemean
    @test res0.av ≈ truemean
end

testEP()

printstyled("All TestEP passed!\n",color=:green,bold=true)

end # end module TestGauge
