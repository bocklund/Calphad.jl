# Implements a `calculate` function for sampling internal degrees of freedom of phases

import IterTools

function sample_phase_constitution(prx::PhaseRecord, pdens)
    phase_constituents = prx.constituent_array
    sublattice_dof = [length(subl) for subl in phase_constituents]
    sampler = halton  # TODO: assumption
    fixed_grid = true  # TODO: assumption

    # Eliminate pure VA endmembers
    vacancy_indices = []
    for subl_idx in 1:length(phase_constituents)
        sublattice = phase_constituents[subl_idx]
        active_in_subl = sort(sublattice)
        subl_va = []
        for sp_idx in 1:length(sublattice)
            if sublattice[sp_idx] == "VA"
                push!(subl_va, sp_idx)
            end # if
        end # for
        push!(vacancy_indices, subl_va)
    end # for

    if length(vacancy_indices) == length(phase_constituents)
        # There is a vacancy in each sublattice, these need to be removed from points
        points = endmember_matrix(sublattice_dof, vacancy_indices)
    else
        points = endmember_matrix(sublattice_dof)
    end # if

    if fixed_grid
        # Sample along endmember edges
        em_pairs = collect(IterTools.subsets(points, 2))
        lingrid1 = LinRange(0, 1, pdens)'
        lingrid2 = reverse!(LinRange(0, 1, pdens))'
        extra_points = [em1*lingrid + em2*lingrid for (em1, em2) in em_pairs]
    end # if

end # function
