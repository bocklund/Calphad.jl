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

end # module
