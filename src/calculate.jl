
# Assume broadcasting
# statevar_dict::Dict
function compute_phase_values(statevar_dict, prx, points)
    # For now, does not produce the
    svd = statevar_dict
    statevars = sort!(collect(keys(svd)))
    statevars_shape = Tuple([length(svd[sv]) for sv in statevars])
    nbpts = size(points)[1]
    vals = fill(NaN, statevars_shape..., nbpts)
    for idx in CartesianIndices(statevars_shape)
        for pt in 1:nbpts
            vals[idx, pt] = prx.obj([svd[statevars[i]][idx[i]] for i in 1:length(idx)]..., points[pt, :]...)
        end # for
    end # for
    return vals
end # function
