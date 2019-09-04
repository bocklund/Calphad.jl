# Precompiling turned off for PyCall. Can be worked around.
# See: https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
__precompile__(false)
module Calphad

# Core structures and costants
include("phase_record.jl")
include("result.jl")
include("core/compositionset.jl")
include("core/constants.jl")

# Files with no internal dependencies
include("core/halton.jl")

# Everything else
include("pycalphad.jl")
include("solver.jl")
include("core/hyperplane.jl")
include("calculate.jl")
include("grids.jl")

export PhaseRecord, Database, Model, calculate

end # module
