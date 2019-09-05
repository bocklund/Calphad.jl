abstract type Result end

struct Grid <: Result
    X::Array{Float64, 2}
    Y::Array{Float64, 2}
    Phase::Array{String, 1}
    output::Array{Float64, 2}
    points::Array{Int64, 1}  # indices
    statevars::Dict{String, Array{Float64, 1}}
end
