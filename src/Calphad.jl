module Calphad
using Symbolics
import LinearAlgebra
using OrderedCollections

export CompSet, PhaseRecord
export R, P, T

include("constants.jl")
include("structs.jl")
include("solver.jl")

end
