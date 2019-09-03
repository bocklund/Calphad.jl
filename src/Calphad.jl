# Precompiling turned off for PyCall. Can be worked around.
# See: https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
__precompile__(false)
module Calphad

include("phase_record.jl")
include("core/compositionset.jl")
include("core/constants.jl")
include("core/halton.jl")
include("pycalphad.jl")
include("solver.jl")
include("core/hyperplane.jl")
include("calculate.jl")
include("grids.jl")

export PhaseRecord, Database, Model

end # module
