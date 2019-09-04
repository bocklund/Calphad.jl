
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

@testset "compute composition" begin
    dbf = Database("data/Cu-Ni.tdb")
    model = Model(dbf, ["CU", "NI"], "LIQUID")
    prx = PhaseRecord(model)
    points = Calphad.sample_phase_constitution(prx, 11)
    X = Calphad.compute_composition(prx, ["CU", "NI"], points)
    # No VA in LIQUID, should be the same as points
    @test all(X .== points)

    # FCC has VA in the second sublattice: (CU, NI):(VA)
    model = Model(dbf, ["CU", "NI", "VA"], "FCC_A1")
    prx = PhaseRecord(model)
    points = Calphad.sample_phase_constitution(prx, 5)
    @test all(size(points) .== (12, 3))
    # Only pure elements have compositions
    X = Calphad.compute_composition(prx, ["CU", "NI"], points)
    @test all(size(X) .== (12, 2))
    @test all(sum(X, dims=2) .≈ ones(12))

    X = Calphad.compute_composition(prx, ["CU", "NI", "CR", "FE"], points)
    @test all(size(X) .== (12, 4))
    # sorted index of CR and FE are 1 and 3. Both columns should be zero.
    @test all(X[:, 1] .== 0.0)
    @test all(X[:, 3] .== 0.0)
    @test all(sum(X, dims=2) .≈ ones(12))

end # begin

@testset "calculate" begin
    dbf = Database("data/Cu-Ni.tdb")
    phases = ["FCC_A1", "LIQUID"]
    prxs = []
    for ph in phases
        model = Model(dbf, ["CU", "NI", "VA"], ph)
        push!(prxs, PhaseRecord(model))
    end # for

    statevars = Dict("T"=>[500.0, 1000.0, 1500.0, 2000.0])
    cr = calculate(prxs, statevars, 11)
    @test isa(cr, Calphad.CalculateResult)
    npts = 48
    @test all(size(cr.X) .== (npts, 2))
    @test all(size(cr.Y) .== (npts, 3))
    @test all(size(cr.Phase) .== (npts,))
    @test all(size(cr.output) .== (4, npts))  # 4 = statevars shape
    @test all(size(cr.statevars["T"]) .== size(statevars["T"]))

    pycalphad_result = [ # pycalphad=0.8
    -17995.39092277  -16427.96788138  -16427.96788138  -16762.79524831  -16886.64602002  -17046.89858755  -17252.83749087  -17500.24908573  -17775.58465488  -18052.93101717  -18287.57204886  -18392.22135906  -17995.39092277  -17862.77142867  -17657.55908095  -16834.42709053  -18338.57811069  -17357.95153724  -17062.86008222  -17137.9987201  -18233.67237248  -16860.94680913  -18277.71775113  -16869.97974943   -9785.71804245   -4405.16813996   -4405.16813996   -5052.08401134   -5397.02601571   -5776.28216043   -6241.93881363   -6803.51045759   -7451.54687413   -8162.12376142   -8892.35211721   -9556.77105334   -9785.71804245   -7666.1799518    -7169.28954376   -5262.58347382   -9785.13608179   -6479.24424949   -5812.73172914   -5983.1088881   -8701.99586306   -5331.90859483   -8855.71772307   -5355.00958336;  -46322.69712689  -44814.72877642  -44814.72877642  -46552.53697683  -47330.63886137  -47839.36526331  -48201.72817811  -48456.8632198  -48606.57251924  -48623.84186158  -48443.73774097  -47916.18298092  -46322.69712689  -48628.15638016  -48557.4315929   -47068.3224136  -47444.16475554  -48327.26479913  -47875.48605727  -48024.30647033  -48516.1655321   -47210.19378732  -48459.42913927  -47254.30332676  -42873.93944689  -37785.75682439  -37785.75682439  -39730.7076689  -40741.3812994   -41512.69305057  -42171.15497874  -42758.1353684  -43277.07310324  -43702.45889957  -43970.8540729   -43933.6761669  -42873.93944689  -43420.66281215  -43067.78116384  -40383.99380627  -43703.80483292  -42439.27077328  -41572.99503272  -41832.39054157  -43921.35467884  -40574.67339671  -43962.7002794   -40635.17536194;  -82059.8329191   -81682.0940662   -81682.0940662   -84546.99006655  -85835.06171426  -86598.07820879  -87035.41728861  -87218.99133359  -87163.93633337  -86841.99206544  -86168.006203    -84929.11078156  -82059.8329191   -87095.44052902  -87217.00773117  -85409.23307158  -84002.00245032  -87146.17619741  -86647.08332191  -86837.98423174  -86387.56994163  -85641.25073355  -86213.27434701  -85712.6122097  -83454.61262034  -79425.24712425  -79425.24712425  -82591.58795726  -84191.34822923  -85278.07060223  -86052.69282071  -86588.43697144  -86901.63103993  -86965.18076067  -86695.09776689  -85879.67803413  -83454.61262034  -86950.4712752   -86797.58284365  -83641.2457765  -85156.0111804   -86317.64539123  -85355.61379158  -85674.63584309  -86806.60856067  -83937.57949591  -86719.32870663  -84030.29566814; -122965.8621675  -124747.36875776 -124747.36875776 -128630.82763497 -130324.3404573  -131239.05120162 -131649.36482305 -131659.47792907 -131297.76556996 -130534.71470019 -129264.9541044  -127212.82502997 -122965.8621675  -131107.67484913 -131495.96061742 -129775.88948144 -125783.35371605 -131700.17528113 -131292.13050558 -131485.97906388 -129658.44592655 -130076.7732596  -129344.98920654 -130168.37796388 -129312.10036086 -127551.56701253 -127551.56701253 -131818.29633189 -133886.14174321 -135167.2732359  -135937.05424256 -136300.56065223 -136287.00955222 -135867.72169523 -134938.15903221 -133223.49597055 -129312.10036086 -136203.67484744 -136339.16380321 -133191.04923308 -131949.89539744 -136145.20857875 -135251.62034307 -135583.86729932 -135242.49661532 -133568.25053099 -135000.93826452 -133684.99241717]
    @test all(pycalphad_result .≈ cr.output)
end # begin
