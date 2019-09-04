
# Assume broadcasting
# statevar_dict::Dict
# TODO: create a compute_phase_values that fills a preallocated view
function compute_phase_values(statevar_dict, prx, points)
    # For now, does not produce the
    svd = statevar_dict
    statevars = sort!(collect(keys(svd)))
    statevars_shape = Tuple([length(svd[sv]) for sv in statevars])
    nbpts = size(points)[1]
    vals = Array{Float64}(undef, statevars_shape..., nbpts)
    for idx in CartesianIndices(statevars_shape)
        for pt in 1:nbpts
            vals[idx, pt] = prx.obj([svd[statevars[i]][idx[i]] for i in 1:length(idx)]..., points[pt, :]...)
        end # for
    end # for
    return vals
end # function

"""

x_i = \\frac{ \\sum_s  (n^s y^s_i) }{ \\sum_s (1 - y^s_{\\mathrm{VA}}) }
"""
function compute_composition(prx, pure_elements::Array{String}, points)
    # Components: points to compute the composition for
    # TODO: Assumes that pure_elements have no VA, and are not species
    pure_elements = sort(pure_elements)
    constituents = prx.constituent_array
    site_ratios = prx.subl_site_ratios
    nelements = length(pure_elements)
    npoints = size(points)[1]
    nsubl = length(site_ratios)
    X = Array{Float64}(undef, npoints, nelements)

    # We only need to create one vacancy normalizing part, so do that first
    va_norm_factor = fill(0.0, npoints)
    dof_offset = 0
    for subl_idx in 1:nsubl
        ns = site_ratios[subl_idx]
        subl_constituents = constituents[subl_idx]
        # get the VA part
        va_idx = findfirst(xx->xx=="VA", subl_constituents)
        if va_idx === nothing
            va_norm_factor .+= ns*1.0
        else
            ys_va = points[:, va_idx+dof_offset]
            va_norm_factor += ns*(1.0 .- ys_va)
        end # if
        dof_offset += length(subl_constituents)
    end # for

    for el_idx in 1:nelements
        el = pure_elements[el_idx]
        comp_amnts = fill(0.0, npoints)
        dof_offset = 0
        for subl_idx in 1:nsubl
            ns = site_ratios[subl_idx]
            subl_constituents = constituents[subl_idx]
            nconstituents = length(subl_constituents)
            for sp_idx in 1:nconstituents
                sp = subl_constituents[sp_idx]
                if sp == el
                    ys_i = points[:, sp_idx+dof_offset]
                    comp_amnts += ns*ys_i
                end # if
            end # for
            dof_offset += length(nconstituents)
        end # for
        X[:, el_idx] = comp_amnts ./ va_norm_factor
    end # for
    return X
end # function

function extend_points(points, max_internal_dof)
    # TODO: eventually we should preallocate the points so this is unnecessary
    # In that case, we can still keep this method around as a no-op
    nbpts = size(points)[1]
    additional_dof = max_internal_dof - size(points)[2]
    if additional_dof > 0
        return [points fill(NaN, nbpts, additional_dof)]
    elseif additional_dof == 0
        return points
    else
        error("Maximum internal degrees of freedom must be greater than or equal to degree of freedom")
    end # if
end # function

function get_nonvacant_pure_elements(species)
    # TODO: assumes species are strings and all pure elements
    elements = unique!(species)
    if "VA" in elements
        deleteat!(elements, findfirst(x->x=="VA", elements))
    end # if
    return sort!(elements)
end # function

function calculate(phase_records, statevars, pdens)
    species = vcat([vcat(prx.constituent_array...) for prx in phase_records]...)
    nonvacant_pure_elements = get_nonvacant_pure_elements(species)
    phase_names = []
    point_nbrs = []
    allpoints = []
    phase_values = []
    max_internal_dof = 0
    compositions = []
    for prx in phase_records
        points = sample_phase_constitution(prx, 11)
        push!(phase_names, prx.name)
        push!(point_nbrs, size(points)[1])
        push!(allpoints, points)
        push!(phase_values, compute_phase_values(statevars, prx, points))
        max_internal_dof = max(size(points)[2], max_internal_dof)
        push!(compositions, compute_composition(prx, nonvacant_pure_elements, points))
    end # for
    site_fracs = vcat([extend_points(pts, max_internal_dof) for pts in allpoints]...)
    composition = vcat(compositions...)
    phases = vcat([fill(phase_names[i], size(allpoints[i])[1]) for i in 1:length(phase_records)]...)
    output = hcat(phase_values...)
    ptsidx = collect(1:size(allpoints)[1])

    return CalculateResult(composition, site_fracs, phases, output, ptsidx, statevars)
end # function
