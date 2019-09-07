mutable struct CompositionSet
	name::String
	dof::Array{Float64}
	NP::AbstractFloat
	energy::AbstractFloat
	phase_record::PhaseRecord
end # struct
phaserecord(x::CompositionSet) = x.phase_record
