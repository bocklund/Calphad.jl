# Precompiling turned off for PyCall. Can be worked around.
# See: https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
module Calphad

using LinearAlgebra
using Symbolics

# Core structures and files no internal dependencies
include("constants.jl")
export R, P, T
include("solver_structs.jl")
export CompSet, PhaseRecord

# Dependent structures
include("symbolic_solver.jl")

# Everything else



end # module
