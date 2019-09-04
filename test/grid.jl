# Tests for sampling grids

@testset "grids" begin
    dbf = Database("data/Cu-Ni.tdb")
    model = Model(dbf, ["CU", "NI"], "LIQUID")
    prx = PhaseRecord(model)
    points = Calphad.sample_phase_constitution(prx, 4)
    @test all(size(points) .== (10, 2))  # (2 endmembers + 4 grid + 4 points, 2 sites)

    pycalphad_result = [ # pycalphad==0.8
    1.00000000e+00 1.00000000e-15;
    1.00000000e-15 1.00000000e+00;
    1.00000000e-15 1.00000000e+00;
    3.33333333e-01 6.66666667e-01;
    6.66666667e-01 3.33333333e-01;
    1.00000000e+00 1.00000000e-15;
    6.30929754e-01 3.69070246e-01;
    5.57885891e-01 4.42114109e-01;
    1.60558422e-01 8.39441578e-01;
    9.46394630e-01 5.36053696e-02
    ]

    @test all(points .≈ pycalphad_result)
end # begin
