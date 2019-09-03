
@testset "pycalphad interface" begin
    # Database loading
    dbf = Database("data/Cu-Ni.tdb")
    # elements:
    @test length(dbf.elements) == 4  # elements: /-, CU, NI, VA
    @test length(dbf._parameters.all()) == 14

    # Model successfully instantiated
    mod = Model(dbf, ["CU", "NI", "VA"], "FCC_A1")
    @test length(mod.pure_elements) == 3  # CU, NI, VA
    @test length(mod.nonvacant_elements) == 2  # CU, NI

    # function can be created
    myfunc, myargs = Calphad.getfunc(mod, "GM")  # TODO: need a helper PhaseRecord that takes a Model instance
    @test all(myargs .== ["T", "FCC_A10CU", "FCC_A10NI", "FCC_A11VA"])
    @test myfunc(500.0, 0.5, 0.5, 1.0) ≈ -17500.249086  # pycalphad verified

    prx = PhaseRecord("FCC_A1", myfunc, myargs, Calphad.active_constituents(mod), Calphad.get_site_ratios(mod))
    @test prx.obj(500.0, 0.5, 0.5, 1.0) ≈ -17500.249086
    @test prx.name == "FCC_A1"
    @test all(prx.constituent_array .== [["CU", "NI"], ["VA"]])
    @test all(prx.args .== ["T", "FCC_A10CU", "FCC_A10NI", "FCC_A11VA"])
    @test all(prx.subl_site_ratios .== [1, 1])
end # begin
