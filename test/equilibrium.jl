import JuMP
@testset "equilibrium" begin
    #dbf = Database("data/Cu-Ni.tdb")
    dbf = Database("data/Cu-Ni.tdb")
    comps = ["CU", "NI"]
    phases = ["A1", "LIQUID"]
    prxs = PhaseRecord[]
    for ph in phases
        model = Model(dbf, comps, ph)
        push!(prxs, PhaseRecord(model, ["T", "P"]))
    end # for
    # Two phase
    conditions = Dict("T"=>[1573.0], "P"=>[101325.0], "X_NI"=>[0.5])  # two phase region
    statevars = Dict([(ky, conditions[ky]) for ky in keys(conditions) if !any(startswith.(ky, ["X_", "MU_"]))])
    grid = calculate(prxs, statevars, 11)
    particular_conds = Dict([(ky, conditions[ky][1]) for ky in keys(conditions)]) # e.g. how these would be found iteratively
    compsets = starting_point(prxs, comps, grid, CartesianIndex(1, 1), particular_conds) # hardcoded index of conditions
    prb = local_equilibrium!(compsets, ["CU", "NI"], particular_conds)
    # pycalphad 0.8 validated
    @test compsets[1].NP ≈ 0.7347958179763151
    @test compsets[2].NP ≈ 0.2652041820236848
    @test all(compsets[1].dof .≈ [0.5375428184637843, 0.46245718153621573])
    @test all(compsets[2].dof .≈ [0.3959808786131095, 0.6040191213868904])
    @test Calphad.energy(prb) ≈ -93526.30696943

    # Single phase
    dbf = Database("data/Cu-Ni.tdb")
    comps = ["CU", "NI"]
    phases = ["A1", "LIQUID"]
    prxs = PhaseRecord[]
    for ph in phases
        model = Model(dbf, comps, ph)
        push!(prxs, PhaseRecord(model, ["T", "P"]))
    end # for
    conditions = Dict("T"=>[1573.0, 1200.0], "P"=>[101325.0], "X_NI"=>[0.5])  # two phase region
    statevars = Dict([(ky, conditions[ky]) for ky in keys(conditions) if !any(startswith.(ky, ["X_", "MU_"]))])
    grid = calculate(prxs, statevars, 11)
    particular_conds = Dict("T"=>1200.0, "P"=>101325.0, "X_NI"=>0.5) # e.g. how these would be found iteratively
    compsets = starting_point(prxs, comps, grid, CartesianIndex(1, 2), particular_conds) # hardcoded index of conditions
    prb = local_equilibrium!(compsets, ["CU", "NI"], particular_conds)
    # NP values are arbitrary because it's the same phase
    @test all(compsets[1].dof .≈ [0.5, 0.5])
    @test all(compsets[2].dof .≈ [0.5, 0.5])
    @test Calphad.energy(prb) ≈ -63185.74156548

end # begin
