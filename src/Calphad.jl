# Precompiling turned off for PyCall. Can be worked around.
# See: https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
__precompile__(false)
module Calphad

# Core structures and files no internal dependencies
include("constants.jl")
include("phase_record.jl")
include("result.jl")
include("sampling.jl")
include("utils.jl")

# Dependent structures
include("compositionset.jl")
include("starting_point.jl")
include("solver.jl")

# Everything else
include("calculate.jl")
include("equilibrium.jl")
include("pycalphad.jl")

export Database, CompositionSet, Grid, Model, PhaseRecord
export calculate, local_equilibrium!, starting_point

end # module
