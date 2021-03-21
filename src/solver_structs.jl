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

function phase_matrix(hess, num_potentials, num_site_fractions, sublattice_dof)
    phase_mat = Matrix{Num}(undef, num_site_fractions + length(sublattice_dof), num_site_fractions + length(sublattice_dof))
    # Fill the upper-left corner with the Hessian terms
    for i in 1:num_site_fractions
        for j in 1:num_site_fractions
            phase_mat[i,j] = hess[num_potentials+i,num_potentials+j]
        end
    end
    # Fill the remaining rows (and the symmetric terms in the columns)
    # using the condition that the sum of site fractions is one
    subl_offset = 0
    for s in 1:length(sublattice_dof)
        # The first columns until the offset are zero
        println("row $(s+num_site_fractions), 1:$subl_offset = 0, $(subl_offset+1):$(subl_offset+sublattice_dof[s]) = 1, $(subl_offset+sublattice_dof[s]+1):end = 0")
        phase_mat[num_site_fractions+s,1:subl_offset] .= 0.0
        # rows in symmetric column
        phase_mat[1:subl_offset,num_site_fractions+s] .= 0.0
        # The columns for the current sublattice are all one
        phase_mat[num_site_fractions+s,subl_offset+1:subl_offset+sublattice_dof[s]] .= 1.0
        # rows in symmetric column
        phase_mat[subl_offset+1:subl_offset+sublattice_dof[s],num_site_fractions+s] .= 1.0
        # Remaining columns are zero
        phase_mat[num_site_fractions+s,subl_offset+sublattice_dof[s]+1:end] .= 0.0
        # rows in symmetric columns,
        # some repeated work here, but it's consistent and logically equal
        phase_mat[subl_offset+sublattice_dof[s]+1:end,num_site_fractions+s] .= 0.0
        subl_offset += sublattice_dof[s]
    end
    println(phase_mat)
    return phase_mat
end


function PhaseRecord(phase_name, obj, mass, potentials, site_fractions, sublattice_dof)
    grad = Symbolics.gradient(obj, [potentials..., site_fractions...])
    hess = Symbolics.hessian(obj, [potentials..., site_fractions...])
    mass_jac = Symbolics.jacobian(mass, [potentials..., site_fractions...])
    pm = phase_matrix(hess, length(potentials), length(site_fractions), sublattice_dof)

    return PhaseRecord(phase_name,
                       obj,
                       grad,
                       hess,
                       mass,
                       mass_jac,
                       potentials,
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

