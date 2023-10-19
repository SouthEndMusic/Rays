module Rays

using Base.Threads
using Combinatorics: combinations
using IterTools: product
using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!

include("camera.jl")
include("shape.jl")
include("intersection.jl")
include("render.jl")

end # Rays
