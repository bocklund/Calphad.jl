struct CompositionSet
	name::String
	dof::Array{Real}
	NP::AbstractFloat
	energy::AbstractFloat
	phase_record::PhaseRecord
end # struct
phaserecord(x::CompositionSet) = x.phase_record
