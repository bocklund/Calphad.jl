

function local_equilibrium!(compsets::Array{CompositionSet, 1}, components::Array{String, 1}, conditions::Dict{String, Float64})
    # TODO: edge cases where compsets of the same phase will raise an error
    # TODO: implicit assumption of P, T conditions here and in solver
    all_phase_records = PhaseRecord[cs.phase_record for cs in compsets]
    problem = build_model(all_phase_records, components)

    # TODO: make these constraints parameters so they may be changed in tight loops without updating the model
    # TODO: warning in filter here for single pair instead of two arguments
    for (cond, val) in statevar_conds(conditions)
        JuMP.@constraint(problem, JuMP.variable_by_name(problem, cond) == val)
    end # for
    for (cond, val) in filter((k, v)->startswith(k, "X_"), conditions)
        # TODO: verify add_composition_constraint handles VA correctly
        add_composition_constraint(problem, all_phase_records, cond[3:end], val, components)
    end # for

    # TODO: support for chempot conditions, how to do this in the problem?
    if any(startswith.(keys(conditions), "MU_"))
        error("Chemical potential conditions not currently supported")
    end # if

    # initial guess
    for i in 1:length(compsets)
        cs = compsets[i]
        phase_name = cs.name
        var = JuMP.variable_by_name(problem, "NP$i")
        JuMP.set_start_value(var, cs.NP)
        carr = cs.phase_record.constituent_array
        JuMP.set_start_value.(internal_dof_variable(problem, cs.phase_record, i), cs.dof)
    end # for
    JuMP.
    JuMP.optimize!(problem)

    # TODO: update compsets; for now, just return the problem
    return problem
end # function
