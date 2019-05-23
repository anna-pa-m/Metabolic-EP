module TestIO
using MetabolicEP, Test, MAT

function testIO()

  S = [1.0 1.0 1.0 -1.0;
       1.0 1.0 1.0 -1.0]
  b=zeros(size(S,1))
  lb=zeros(size(S,2))
  ub=ones(size(S,2));
  rxnNames = ["A","B","C","D"]
  testdict = Dict{String,Any}("test" => Dict{String,Any}("S"=>S,"b"=>b,"rxnNames"=>rxnNames, "lb"=>lb,"ub"=>ub));
  matwrite("test.mat", testdict)
  X = ReadMatrix("test.mat")
  @test rxnNames == X.rxnNames
  @test b == X.b
  @test lb == X.lb
  @test ub == X.ub
  @test S == X.S
end

testIO()
rm("./test.mat",force=true)
printstyled("All TestIO passed!\n",color=:green,bold=true)

end # end module TestGauge
