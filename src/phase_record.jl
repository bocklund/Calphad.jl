#=
PhaseRecords should have a 1:1 correspondence with each phase.

They should contain the Gibbs energy function (and any derivatives)
as well as information about how the internal constraints for that phase
are constructed and how the parts of that phase can be used externally.

For example:
1. The internal degrees of freedom need to be set up so that they can be
   splatted into the objective function and built into constraints
2. The sublattice constituents and site ratios must be accounted so the
   number of moles can be counted.
=#


# AbstractPhaseRecord should have one method in the interface,
# `phaserecord`, which retrives the phase record

abstract type AbstractPhaseRecord end

struct PhaseRecord <: AbstractPhaseRecord
	name::String
	obj::Function
#	jac::Function
#   hess ::Function
	args::Array{String, 1}
	constituent_array::Array{Array{String, 1}}
	subl_site_ratios::Array{AbstractFloat}
end
phaserecord(x::PhaseRecord) = x
