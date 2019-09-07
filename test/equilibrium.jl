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
    println(JuMP.value.(JuMP.variable_by_name.(prb, ["NP1", "NP2", "INTERNAL_DOF1[1]", "INTERNAL_DOF1[2]", "INTERNAL_DOF2[1]", "INTERNAL_DOF2[2]"])))


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
    println(JuMP.value.(JuMP.variable_by_name.(prb, ["NP1", "NP2", "INTERNAL_DOF1[1]", "INTERNAL_DOF1[2]", "INTERNAL_DOF2[1]", "INTERNAL_DOF2[2]"])))


end # begin