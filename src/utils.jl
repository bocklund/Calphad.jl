function get_nonvacant_pure_elements(species)
    # TODO: assumes species are strings and all pure elements
    elements = unique!(species)
    if "VA" in elements
        deleteat!(elements, findfirst(x->x=="VA", elements))
    end # if
    return sort!(elements)
end # function

statevarfilter(x) = !any(startswith.(x, ["X_", "MU_"]))
statevarfilter(kv::Pair) = statevarfilter(kv[1])
statevar_conds(conditions) = filter(statevarfilter, conditions)
non_statevar_conds(conditions) = filter(xx->!statevarfilter(xx), conditions)
