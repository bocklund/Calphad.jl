struct PhaseRecord
    phase_name
    obj
    grad
    hess
    mass
    mass_jac
    state_variables
    site_fractions
    phase_matrix
    inv_phase_matrix
end

function phase_matrix(hess, num_state_variables, num_site_fractions)
    num_internal_constraints = 1  # TODO: hardcoded for now
    num_internal_dof = num_site_fractions
    phase_mat = Matrix{Num}(undef, num_internal_dof + num_internal_constraints, num_internal_dof + num_internal_constraints)
    for i in 1:num_internal_dof
        for j in 1:num_internal_dof
            phase_mat[i,j] = hess[num_state_variables+i,num_state_variables+j]
        end
    end
    for i in 1:num_internal_constraints
        for j in 1:num_internal_constraints
            phase_mat[num_internal_dof+i,:] .= 1.0
            phase_mat[:,num_internal_dof+j] .= 1.0
            phase_mat[num_internal_dof+i, num_internal_dof+i] = 0.0
        end
    end
    return phase_mat
end


function PhaseRecord(phase_name, obj, mass, state_variables, site_fractions)
    grad = Symbolics.gradient(obj, [state_variables..., site_fractions...])
    hess = Symbolics.hessian(obj, [state_variables..., site_fractions...])
    mass_jac = Symbolics.jacobian(mass, [state_variables..., site_fractions...])
    pm = phase_matrix(hess, length(state_variables), length(site_fractions))

    return PhaseRecord(phase_name,
                       obj,
                       grad,
                       hess,
                       mass,
                       mass_jac,
                       state_variables,
                       site_fractions,
                       pm,
                       LinearAlgebra.inv(pm),
                       )
end

mutable struct CompSet
    phase_rec::PhaseRecord
    Y
    ℵ
end

