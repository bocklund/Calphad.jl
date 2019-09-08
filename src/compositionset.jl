mutable struct CompositionSet
	name::String
	dof::Array{Float64}
	NP::AbstractFloat
	energy::AbstractFloat
	phase_record::PhaseRecord
end # struct
phaserecord(x::CompositionSet) = x.phase_record

function unique_compsets(compsets::Array{CompositionSet, 1})
	newcs = [compsets[1]]
	if length(compsets) == 1
		return newcs
	end # for
	for proposed in compsets[2:end]
		matched = false
		for cs in newcs
			if (proposed.name == cs.name) && (proposed.dof ≈ cs.dof)
				matched = true
				break
			end # if
		end # for
		if !matched
			push!(newcs, proposed)
		end # if
	end # for
	return newcs
end # for
