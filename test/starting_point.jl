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
end # begin
