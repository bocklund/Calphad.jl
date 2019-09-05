function get_nonvacant_pure_elements(species)
    # TODO: assumes species are strings and all pure elements
    elements = unique!(species)
    if "VA" in elements
        deleteat!(elements, findfirst(x->x=="VA", elements))
    end # if
    return sort!(elements)
end # function
