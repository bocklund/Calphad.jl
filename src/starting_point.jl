import LinearAlgebra

# Solve Ax=b
function solve(A, b)
    x = try
        A\b
    catch sing
        if isa(sing, LinearAlgebra.SingularException)
            return -1.0e19
        else
            rethrow()
        end # if
    end
    return x
end # function

"""

    hyperplane(compositions, energies, target_composition, chemical_potentials, total_moles,
                    fixed_chempot_indices, fixed_comp_indices, result_fractions, result_simplex)

Find chemical potentials which approximate the tangent hyperplane at the given composition.

Energies and chemical potentials should come from one set of state variables.

# Arguments
- `compositions::Array{AbstractFloat,2}`: (M, N) normalized compositions corresponding to
  the energy surface
- `energies::Array{AbstractFloat,1}`: (M,) array of molar Gibbs energy surface
- `target_composition::Array{AbstractFloat,1}`: (N,) array of the desired composition
  for the hyperplane.
- `chemical_potentials::Array{AbstractFloat,1}`: (N,) array of chemical potentials
  (including any that are fixed), will be updated with the found hyperplane.
- `total_moles::AbstractFloat`: Total number of moles in the system.
- `fixed_chempot_indices::Array{Integer,1}`: Indices of fixed chemical potentials, variable
  shape from (0,) to (N-1,)
- `fixed_comp_indices::Array{Integer,1}`: Indices of compositions fixed by the conditions,
  variable shape from (0,) to (N-1,)
- `result_fractions::Array{Float64,1}`: Weighted amounts of the compositions making up the
  hyperplane simplex, shape of (P,), will be overwritten. Output sums to 1.0.
- `result_simplex::Array{Integer,1}`: Indices of the points making up the hyperplane
  simplex. Shape of (P,). Will be overwritten.

# Examples
None yet.

# Notes
M: number of energy points that have been sampled
N: number of components
P: N+1, max phases by Gibbs phase rule that we can find in a point calculation
"""
function hyperplane(compositions, energies, target_composition, chemical_potentials, total_moles,
                    fixed_chempot_indices, fixed_comp_indices, result_fractions, result_simplex)
    num_components = size(compositions)[2]
    num_points = length(energies)
    num_fixed_chempots = size(fixed_chempot_indices)[1]
    simplex_size = num_components - num_fixed_chempots
    # composition index of -1 indicates total number of moles, i.e., N=1 condition
    included_composition_indices = [sort(fixed_comp_indices)..., -1]
    best_guess_simplex = sort(setdiff!(collect(1:num_components), fixed_chempot_indices))
    free_chempot_indices = best_guess_simplex[:]
    candidate_simplex = best_guess_simplex[:]
    trial_simplices = Array{Int}(undef, simplex_size, simplex_size)
    fractions = Array{Float64}(undef, simplex_size, simplex_size)
    driving_forces = Array{Float64}(undef, num_points)
    for i in 1:simplex_size
        trial_simplices[i, :] = best_guess_simplex
    end
    trial_matrix = Array{Float64}(undef, simplex_size, simplex_size, simplex_size)
    candidate_tieline = Array{Float64}(undef, simplex_size, simplex_size)
    candidate_energies = Array{Float64}(undef, simplex_size)
    candidate_potentials = Array{Float64}(undef, simplex_size)
    smallest_fractions = Array{Float64}(undef, simplex_size)
    contig_trial = Array{Float64}(undef, simplex_size, simplex_size)
    saved_trial = 0

    max_iterations = 1000
    for _ in 1:max_iterations
        for trial_idx in 1:simplex_size
            for comp_idx in 1:simplex_size
                ici = included_composition_indices[comp_idx]
                for simplex_idx in 1:simplex_size
                    if ici > 0
                        trial_matrix[comp_idx, simplex_idx, trial_idx] = compositions[trial_simplices[trial_idx, simplex_idx], ici]
                    else
                        # ici = -1, refers to N=1 condition
                        trial_matrix[comp_idx, simplex_idx, trial_idx] = 1 # 1 mole-formula per formula unit of a phase
                    end # if
                end # for
            end # for
        end # for

        for trial_idx in 1:simplex_size
            contig_trial[:] = trial_matrix[:, :, trial_idx]
            for simplex_idx in 1:simplex_size
                ici = included_composition_indices[simplex_idx]
                if ici > 0
                    fractions[trial_idx, simplex_idx] = target_composition[ici]
                else
                    # ici = -1, refers to N=1 condition
                    fractions[trial_idx, simplex_idx] = total_moles
                end # if
            end # for
            fractions[trial_idx, :] = solve(contig_trial, fractions[trial_idx, :])
            smallest_fractions[trial_idx] = minimum(fractions[trial_idx, :])
        end # for
        # Choose simplex with the largest smallest-fraction
        saved_trial = argmax(smallest_fractions)
        if smallest_fractions[saved_trial] < -simplex_size
            break
        end # if
        # Should be exactly one candidate simplex
        candidate_simplex = trial_simplices[saved_trial, :]
        for i in 1:size(candidate_simplex)[1]
            idx = candidate_simplex[i]
            ici = 0
            for ici in 1:length(free_chempot_indices)
                chempot_idx = free_chempot_indices[ici]
                candidate_tieline[i, ici] = compositions[idx, chempot_idx]
            end # for
            candidate_potentials[i] = energies[idx]
            for j in fixed_chempot_indices
                candidate_potentials[i] -= chemical_potentials[j] * compositions[idx, j]
            end # for
        end # for
        candidate_potentials = solve(candidate_tieline, candidate_potentials)
        if candidate_potentials[1] == -1.0e19
            break
        end # if

        # Compute the driving forces by subtracting out the fixed (caller-provided) and
        # free (candidate) chemical potentials from the energy surface. A negative driving
        # force indicates that a more stable composition set exists.
        driving_forces[:] = energies
        for ici in 1:length(free_chempot_indices)
            chempot_idx = free_chempot_indices[ici]
            for idx in 1:size(driving_forces)[1]
                driving_forces[idx] -= candidate_potentials[ici] * compositions[idx, chempot_idx]
            end # for
        end #
        for chempot_idx in fixed_chempot_indices
            for idx in 1:size(driving_forces)[1]
                driving_forces[idx] -= chemical_potentials[chempot_idx] * compositions[idx, chempot_idx]
            end # for
        end # for
        best_guess_simplex[:] = candidate_simplex
        for i in 1:size(trial_simplices)[1]
            trial_simplices[i, :] = best_guess_simplex
        end # for
        # Trial simplices will be the current simplex with each vertex replaced by the
        # trial point. Exactly one of those simplices will contain a given test point,
        # excepting edge cases.
        min_df = argmin(driving_forces)
        for i in 1:simplex_size
            trial_simplices[i, i] = min_df
        end # for
        if driving_forces[min_df] > -1e-8
            break
        end # if
    end # for
    out_energy = 0
    for i in 1:size(best_guess_simplex)[1]
        idx = best_guess_simplex[i]
        out_energy = out_energy + fractions[saved_trial, i] * energies[idx]
    end # for
    result_fractions[1:simplex_size] = fractions[saved_trial, :]
    for ici in 1:length(free_chempot_indices)
        chempot_idx = free_chempot_indices[ici]
        chemical_potentials[chempot_idx] = candidate_potentials[ici]
    end # for
    result_simplex[1:simplex_size] = best_guess_simplex
    out_energy
end # function


"""
    starting_point(phase_records, grid, target_composition, target_chempots, fixed_compositions_idx, fixed_chempots_idx)
"""
function starting_point(phase_records::Array{PhaseRecord, 1}, grid::Grid, statevar_index::CartesianIndex,
    target_composition::Array{Float64, 1}, target_chempots::Array{Float64, 1},
    fixed_compositions_idx::Array{Int, 1}, fixed_chempots_idx::Array{Int, 1})
    prx_dict = Dict([(prx.name, prx) for prx in phase_records])
    npts = length(grid.points)
    ncomps = length(target_composition)
    simplex_size = ncomps
    total_moles = 1.0;  # TODO: assumption
    result_fractions = zeros(Float64, simplex_size)
    # result_simplex = zeros(Int, simplex_size)
    # filled with larger than points so we know these are the invalid ones
    result_simplex = fill(npts+1, simplex_size)
    hyperplane(grid.X, grid.output[statevar_index, :], target_composition, target_chempots, total_moles, fixed_chempots_idx, fixed_compositions_idx, result_fractions, result_simplex)
    compsets = CompositionSet[]
    for i in 1:simplex_size
        pts_idx = result_simplex[i]
        if pts_idx <= npts  # invalid points are npts+1
            phase_name = grid.Phase[pts_idx]
            cs = CompositionSet(phase_name, grid.Y[pts_idx, :], result_fractions[i], grid.output[pts_idx], prx_dict[phase_name])
            push!(compsets, cs)
        end # if
    end # for
    return compsets
end # function

#comps = get_nonvacant_pure_elements(vcat([vcat(prx.constituent_array...) for prx in phase_records]...)
# TODO: starting_point, but with conditions Dict instead of target_(comps/chempots) and fixed_(comps/chempots)
function starting_point(phase_records::Array{PhaseRecord, 1}, components::Array{String, 1}, grid::Grid, statevar_index::CartesianIndex, conditions::Dict{String, Float64})
    components = sort(components)
    ncomps = length(components)
    indep_comp_conds = [c for c in keys(conditions) if any(startswith.(c, ["X_", "MU_"]))]
    fixed_comps = sort!([c[3:end] for c in indep_comp_conds if startswith(c, "X_")])
    fixed_chempots = sort!([c[4:end] for c in indep_comp_conds if startswith(c, "MU_")])
    fixed_comp_indices = Int[findfirst(xx->xx==c, components) for c in fixed_comps]
    fixed_chempot_indices = Int[findfirst(xx->xx==c, components) for c in fixed_chempots]

    # The target composition and chempots must match the total number of components for hyperplane
    # Replace the fixed values with those in the conditions, leaving the rest undefined
    target_comp = Array{Float64}(undef, ncomps)
    for i in 1:length(fixed_comp_indices)
        comp_idx = fixed_comp_indices[i]
        comp = components[comp_idx]
        target_comp[comp_idx] = conditions["X_$comp"]
    end # for
    target_chempots = Array{Float64}(undef, ncomps)
    for i in 1:length(fixed_chempot_indices)
        comp_idx = fixed_chempot_indices[i]
        comp = components[comp_idx]
        target_comp[comp_idx] = conditions["MU_$comp"]
    end # for
    # phase_records, grid, target_composition, target_chempots, fixed_compositions, fixed_chempots
    return starting_point(phase_records, grid, statevar_index, target_comp, target_chempots, fixed_comp_indices, fixed_chempot_indices)
end # function
