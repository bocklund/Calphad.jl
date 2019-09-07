# Construction grids of internal degrees of freedom

import IterTools


# First 50 multipliers from
# Chi, H., Mascagni, M., & Warnock, T. (2005). On the optimal Halton sequence.
# Mathematics and Computers in Simulation, 70(1), 9-21. doi:10.1016/j.matcom.2005.03.004
# The rest are generated from the procedure discussed in the paper.
# Note that the choice of multipliers is not unique even when using their procedure.
# Some of their multpliers (e.g. for primes 11 and 29) are not primitive roots.
# It's unclear why they work so well despite that, but they do.
    const _scrambling_multipliers = [(2, 1), (3, 2), (5, 2), (7, 5), (11, 3), (13, 7), (17, 3),
                               (19, 10), (23, 18), (29, 11), (31, 17), (37, 5), (41, 17),
                               (43, 26), (47, 40), (53, 14), (59, 40), (61, 44), (67, 12),
                               (71, 31), (73, 45), (79, 70), (83, 8), (89, 38), (97, 82),
                               (101, 8), (103, 12), (107, 38), (109, 47), (113, 70),
                               (127, 29), (131, 57), (137, 97), (139, 110), (149, 32),
                               (151, 48), (157, 84), (163, 124), (167, 155), (173, 26),
                               (179, 69), (181, 83), (191, 157), (193, 171), (197, 8),
                               (199, 32), (211, 112), (223, 205), (227, 15), (229, 31),
                               (233, 61), (239, 112), (241, 127), (251, 212), (257, 7),
                               (263, 57), (269, 108), (271, 120), (277, 178), (281, 210),
                               (283, 234), (293, 34), (307, 161), (311, 199), (313, 219),
                               (317, 255), (331, 63), (337, 120), (347, 218), (349, 237),
                               (353, 278), (359, 341), (367, 58), (373, 118), (379, 176),
                               (383, 218), (389, 282), (397, 369), (401, 12), (409, 93),
                               (419, 194), (421, 221), (431, 325), (433, 350), (439, 419),
                               (443, 21), (449, 86), (457, 173), (461, 219), (463, 241),
                               (467, 288), (479, 425), (487, 28), (491, 78), (499, 168),
                               (503, 215), (509, 283), (521, 429), (523, 455), (541, 138),
                               (547, 211), (557, 335), (563, 410), (569, 485), (571, 510),
                               (577, 13), (587, 134), (593, 210), (599, 286), (601, 309),
                               (607, 383), (613, 463), (617, 521), (619, 539), (631, 75),
                               (641, 205), (643, 232), (647, 282), (653, 363), (659, 442),
                               (661, 467), (673, 635), (677, 17), (683, 92), (691, 202),
                               (701, 334), (709, 444), (719, 584), (727, 701), (733, 52),
                               (739, 136), (743, 191), (751, 310), (757, 388), (761, 446),
                               (769, 560), (773, 615), (787, 45), (797, 183), (809, 362),
                               (811, 389), (821, 538), (823, 567), (827, 624), (829, 657),
                               (839, 811), (853, 173), (857, 236), (859, 264), (863, 325),
                               (877, 538), (881, 600), (883, 632), (887, 695), (907, 105),
                               (911, 168), (919, 292), (929, 444), (937, 569), (941, 636),
                               (947, 735), (953, 828), (967, 89), (971, 155), (977, 251),
                               (983, 347), (991, 474), (997, 573), (1009, 772), (1013, 838),
                               (1019, 940), (1021, 975), (1031, 112), (1033, 146),
                               (1039, 243), (1049, 407), (1051, 442), (1061, 608),
                               (1063, 643), (1069, 747), (1087, 1055), (1091, 37),
                               (1093, 65), (1097, 138), (1103, 234), (1109, 334),
                               (1117, 474), (1123, 574), (1129, 674), (1151, 1067),
                               (1153, 1102), (1163, 117), (1171, 255), (1181, 432),
                               (1187, 540), (1193, 645), (1201, 786), (1213, 1002),
                               (1217, 1076), (1223, 1189), (1229, 71), (1231, 105),
                               (1237, 212), (1249, 430), (1259, 609), (1277, 939),
                               (1279, 977), (1283, 1048), (1289, 1166), (1291, 1202),
                               (1297, 17), (1301, 90), (1303, 126), (1307, 200),
                               (1319, 421), (1321, 456), (1327, 568), (1361, 1213),
                               (1367, 1331), (1373, 75), (1381, 224), (1399, 565),
                               (1409, 756), (1423, 1030), (1427, 1108), (1429, 1148),
                               (1433, 1226), (1439, 1345), (1447, 57), (1451, 134),
                               (1453, 172), (1459, 287), (1471, 520), (1481, 715),
                               (1483, 756), (1487, 836), (1489, 876), (1493, 954),
                               (1499, 1073), (1511, 1316), (1523, 41), (1531, 195),
                               (1543, 433), (1549, 554), (1553, 633), (1559, 755),
                               (1567, 919), (1571, 1000), (1579, 1164), (1583, 1245),
                               (1597, 1541), (1601, 23), (1607, 141), (1609, 183),
                               (1613, 262), (1619, 383), (1621, 425), (1627, 546),
                               (1637, 753), (1657, 1172), (1663, 1295), (1667, 1383),
                               (1669, 1427), (1693, 247), (1697, 331), (1699, 371),
                               (1709, 582), (1721, 835), (1723, 880), (1733, 1091),
                               (1741, 1261), (1747, 1392), (1753, 1527), (1759, 1655),
                               (1777, 275), (1783, 402), (1787, 486), (1789, 534),
                               (1801, 790), (1811, 1007), (1823, 1270), (1831, 1446),
                               (1847, 1805), (1861, 260), (1867, 392), (1871, 476),
                               (1873, 518), (1877, 605), (1879, 654), (1889, 874),
                               (1901, 1141), (1907, 1276), (1913, 1412), (1931, 1822),
                               (1933, 1868), (1949, 288), (1951, 331), (1973, 827),
                               (1979, 963), (1987, 1143)]
const scrambler = Dict([(Float64(s[1]), Float64(s[2])) for s in _scrambling_multipliers])
const _primes = Float64.([s[1] for s in _scrambling_multipliers])

# TODO: overloaded methods where scrambling is not true and primes are passed
# Maybe: YAGNI
function halton(dim, nbpts)

    # dim must be less than or equal to number of scrambling multiplers
    result = Array{Float64}(undef, nbpts, dim);
    dim_primes = _primes[1:dim]
    for i in 1:dim
        num_powers = ceil(Int64, log(nbpts+1)/log(dim_primes[i]))
        powers = dim_primes[i].^(-(1:num_powers))
        radix_vector = dim_primes[i].^(-(0:num_powers-1))
        sum_matrix = Float64.(1:nbpts) * radix_vector'
        mod_matrix = mod.(Float64(scrambler[dim_primes[i]]) * floor.(sum_matrix), dim_primes[i])
        result[:,i] = mod_matrix * powers
    end
    return result
end

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