
# Assume broadcasting
# statevar_dict::Dict
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

function extend_points(points, max_internal_dof)
    # TODO: eventually we should preallocate the points so this is unnecessary
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
