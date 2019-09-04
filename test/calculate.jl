
@testset "compute phase values" begin
    # Compute phase values produces correct output
    dbf = Database("data/Cu-Ni.tdb")
    model = Model(dbf, ["CU", "NI"], "LIQUID")
    prx = PhaseRecord(model)
    points = Calphad.sample_phase_constitution(prx, 11)
    statevars = Dict("T"=>[300.0, 400.0, 500.0])
    phase_vals = Calphad.compute_phase_values(statevars, prx, points)
    pycalphad_values = [ # pycalphad==0.8
        1.66677137e+02  5.52772990e+03  5.52772990e+03  5.41367242e+03  5.34866747e+03  5.13987797e+03  4.76498784e+03  4.22722410e+03  3.54119261e+03  2.73018369e+03  1.82886647e+03  8.97453204e+02  1.66677137e+02  3.30236237e+03  3.84688988e+03  5.38855046e+03  4.92123902e+02  4.54605455e+03  5.11414630e+03  4.98460884e+03  2.07197776e+03  5.37056212e+03  1.87619809e+03  5.36334487e+03;
        -4.48860942e+03  8.88848483e+02  8.88848483e+02  5.07696152e+02  3.02057018e+02  7.36854007e+00 -4.13570513e+02 -9.63903862e+02 -1.63160347e+03 -2.39306203e+03 -3.20950048e+03 -4.00808224e+03 -4.48860942e+03 -1.85854102e+03 -1.33734584e+03  3.89482331e+02 -4.32523823e+03 -6.41984932e+02 -2.37794976e+01 -1.73992075e+02 -2.99259168e+03  3.45689241e+02 -3.16748356e+03  3.30485078e+02;
        -9.78571804e+03 -4.40516814e+03 -4.40516814e+03 -5.05208401e+03 -5.39702602e+03 -5.77628216e+03 -6.24193881e+03 -6.80351046e+03 -7.45154687e+03 -8.16212376e+03 -8.89235212e+03 -9.55677105e+03 -9.78571804e+03 -7.66617995e+03 -7.16928954e+03 -5.26258347e+03 -9.78513608e+03 -6.47924425e+03 -5.81273173e+03 -5.98310889e+03 -8.70199586e+03 -5.33190859e+03 -8.85571772e+03 -5.35500958e+03
    ]
    @test all(phase_vals .≈ pycalphad_values)
end # begin

@testset "extend points" begin
    dbf = Database("data/Cu-Ni.tdb")
    model = Model(dbf, ["CU", "NI"], "LIQUID")
    prx = PhaseRecord(model)
    points = Calphad.sample_phase_constitution(prx, 11)
    @test all(size(points) .== (24, 2))
    extended_points = Calphad.extend_points(points, 5)
    @test all(size(extended_points) .== (24, 5))
end # begin
