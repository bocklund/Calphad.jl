# Precompiling turned off for PyCall. Can be worked around.
# See: https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
__precompile__(false)
module Calphad

# Core structures and files no internal dependencies
include("constants.jl")
include("hyperplane.jl")
include("phase_record.jl")
include("result.jl")
include("sampling.jl")
include("utils.jl")

# Dependent structures
include("compositionset.jl")

# Everything else
include("calculate.jl")
include("pycalphad.jl")
include("solver.jl")

export Database, Grid, Model, PhaseRecord, calculate

end # module
