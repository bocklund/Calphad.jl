# Construction grids of internal degrees of freedom

import IterTools

"""
Accept the number of components in each sublattice.
Return a matrix corresponding to the compositions of all endmembers.

# Arguments
-`dof::Array(Int,1)`: Number of components in each sublattice.
-`vacancy_indices::Array(Int,1)`: Index of vacancy in each sublattice. Assumes that each
  sublattice can only have one vacancy.

# Examples
For a sublattice configuration like: `(AL, NI, VA):(AL, NI, VA):(VA)`
```
endmember_matrix([3,3,1], vacancy_indices=[3, 3, 1])
```
"""
function endmember_matrix(dof, vacancy_indices=[])
    nsubl = length(dof)
    subl_offset = cumsum(vcat([0], dof))[1:end-1]
    # Construct the endmember matrix by a Cartesian product of rows in identity matrices
    # There is one matrix for each sublattice, where each row is one pure species
    eye(n) = [[i == j ? 1.0 : 0.0 for j in 1:n] for i in 1:n]
    res_matrix = permutedims(hcat([vcat(em...) for em in collect(Base.product([eye(n) for n in dof]...))[:]]...))
    # Remove the one vacancy only endmember, if it exists. We essentially copy the array
    # here to remove just one endmember, but it's likely rare that we hit this code path.
    if length(vacancy_indices) == nsubl
        purevaidx = vacancy_indices + subl_offset
        res_matrix = permutedims(hcat([res_matrix[em, :] for em in 1:size(res_matrix)[1] if sum(res_matrix[em, purevaidx] .== 1) != nsubl]...))
    end
    # Normalize site fractions to the numerical limit
    res_matrix[res_matrix .≈ 0] .= MIN_SITE_FRACTION
    for sublidx in 1:length(dof)
        nspecies = dof[sublidx]
        subl = subl_offset[sublidx] .+ (1:nspecies)
        @views res_matrix[:, subl][res_matrix[:, subl] .≈ 1] .= 1.0-MIN_SITE_FRACTION*float(nspecies-1)
    end
    return res_matrix
end

"""
    grid_sample(endmembers, density)

Sample between pairs of endmembers (# endmembers, # internal dof) on a regular grid.
These constitution space edges are often the equilibrium points
"""
function grid_sample(endmembers, density)
    ndof = size(endmembers,2)
    nendmembers = size(endmembers,1)
    em_pairs = collect(IterTools.subsets([endmembers[i,:] for i in 1:nendmembers], 2))
    npairs = size(em_pairs,1)
    line = LinRange(0.0, 1.0, density)
    npts = length(line)  # points per endmember pair
    pts = Array{Float64}(undef, npts*npairs, ndof)
    offset = 0
    for (em1, em2) in em_pairs
        idx = offset.+(1:npts)
        pts[idx, :] = line*reshape(em1, 1, ndof) + reverse(line)*reshape(em2, 1, ndof)
        offset += npts
    end # for
    return pts
end # function

"""
    point_sample(dof, density)

Sample each degree of freedom from a scrambled Halton sequence, normalizing appropriately.

# Arguments
- `dof::Array{Integer,1}`: Number of components in each sublattice
- `density::Integer`: Number of points to sample on each degree of freedom

# Examples
```
point_sample([2, 1], 10)
```
"""
function point_sample(dof::Array{Int,1}, density::Integer)
    nconstituents = sum(dof)
    nsublattices = length(dof)
    ndof = nconstituents - nsublattices
    if ndof == 0
        return ones(1, nsublattices)
    end # if
    pts = halton(nconstituents, density*ndof)
    # Convert low-discrepancy sequence to normalized exponential to distribute the points
    # uniformly over the simplices. For details, see the beginning of section 4 in
    # Otis, Emelianenko, Liu, Comput. Mater. Sci. 130 (2017) 282–291. doi:10.1016/j.commatsci.2017.01.019.
    pts = -log.(pts)
    start_idx = 1 # starting index of next sublattice
    for nsublconsts in dof
        end_idx = start_idx + (nsublconsts-1)
        cur_subl = start_idx:end_idx
        pts[:, cur_subl] ./= sum(pts[:, cur_subl], dims=2)
        start_idx = end_idx + 1
    end # for
    return pts
end # function


function sample_phase_constitution(prx::PhaseRecord, pdens)
    phase_constituents = prx.constituent_array
    sublattice_dof = [length(subl) for subl in phase_constituents]
    sampler = halton  # TODO: assumption
    fixed_grid = true  # TODO: assumption

    # Eliminate pure VA endmembers
    va_idx = []
    for subl_idx in 1:length(phase_constituents)
        sublattice = phase_constituents[subl_idx]
        active_in_subl = sort(sublattice)
        for sp_idx in 1:length(sublattice)
            if sublattice[sp_idx] == "VA"
                push!(va_idx, sp_idx)
            end # if
        end # for
    end # for

    # Endmember points
    if length(va_idx) == length(phase_constituents)
        # There is a vacancy in each sublattice, these need to be removed from points
        points = endmember_matrix(sublattice_dof, va_idx)
    else
        points = endmember_matrix(sublattice_dof)
    end # if

    # Sample along endmember edges
    if fixed_grid
        points = [points; grid_sample(points, pdens)]
    end # if

    # Randomly sample the space
    points = [points; point_sample(sublattice_dof, pdens)]

    return points
end # function
