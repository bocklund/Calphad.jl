module Calphad
using Symbolics
import LinearAlgebra

export CompSet, PhaseRecord
export R, P, T

include("constants.jl")
include("structs.jl")
include("solver.jl")

end
