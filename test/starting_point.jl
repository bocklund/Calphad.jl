@testset "starting point" begin
    dbf = Database("data/Cu-Ni.tdb")
    comps = ["CU", "NI", "VA"]
    phases = ["FCC_A1", "LIQUID"]
    prxs = PhaseRecord[]
    for ph in phases
        model = Model(dbf, comps, ph)
        push!(prxs, PhaseRecord(model))
    end # for
    conditions = Dict("T"=>[500.0], "X_NI"=>[0.5])
    statevars = Dict([(ky, conditions[ky]) for ky in keys(conditions) if !any(startswith.(ky, ["X_", "MU_"]))])
    grid = calculate(prxs, statevars, 11)
    particular_conds = Dict([(ky, conditions[ky][1]) for ky in keys(conditions)]) # e.g. how these would be found iteratively
    compsets = starting_point(prxs, comps, grid, CartesianIndex(1), particular_conds)
    @test all(typeof.(compsets) .== CompositionSet)
    @test length(compsets) == 2
    cs1 = compsets[1]
    cs2 = compsets[2]
	@test cs1.name == "FCC_A1"
	@test cs2.name == "FCC_A1"
	# all values from pycalphad 0.8
	@test all(cs1.dof .≈ [0.7737056144690831, 0.22629438553091683, 1.0])
	@test all(cs2.dof .≈ [0.1, 0.9, 1.0])
	@test cs1.NP ≈ 0.5937311362845352
	@test cs2.NP ≈ 0.4062688637154648
	@test cs1.energy ≈ -18233.67237248
	@test cs2.energy ≈ -16762.79524831
	@test cs1.phase_record === prxs[1]
	@test cs2.phase_record === prxs[1]
end # begin
