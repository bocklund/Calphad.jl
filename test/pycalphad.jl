
@testset "pycalphad PhaseRecord interface" begin
    # Database loading
    dbf = Database("data/Cu-Ni.tdb")
    # elements:
    @test length(dbf.elements) == 4  # elements: /-, CU, NI, VA
    @test length(dbf._parameters.all()) == 14

    # Model successfully instantiated
    model = Model(dbf, ["CU", "NI", "VA"], "FCC_A1")
    @test length(model.pure_elements) == 3  # CU, NI, VA
    @test length(model.nonvacant_elements) == 2  # CU, NI

    prx = PhaseRecord(model)
    @test prx.obj(500.0, 0.5, 0.5, 1.0) ≈ -17500.249086
    @test prx.name == "FCC_A1"
    @test all(prx.constituent_array .== [["CU", "NI"], ["VA"]])
    @test all(prx.args .== ["T", "FCC_A10CU", "FCC_A10NI", "FCC_A11VA"])
    @test all(prx.subl_site_ratios .== [1, 1])
end # begin
