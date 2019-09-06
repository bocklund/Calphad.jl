
@testset "state variable utils" begin
    conds = Dict("T"=> [300.0, 400.0], "P"=>101325.0, "X_NI"=>0.5, "MU_CU"=> 100000.0)
    svs = Calphad.statevar_conds(conds)
    @test "T" in keys(svs)
    @test "P" in keys(svs)
    @test length(svs) == 2

    nsvs = Calphad.non_statevar_conds(conds)
    @test "X_NI" in keys(nsvs)
    @test "MU_CU" in keys(nsvs)
    @test length(nsvs) == 2
end # begin
