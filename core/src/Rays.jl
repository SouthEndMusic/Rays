module Rays

import SimpleDirectMediaLayer as SDL

using Base.Threads
using Combinatorics: combinations
using IterTools: product
using SimpleDirectMediaLayer.LibSDL2
using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!

include("utils.jl")
include("camera.jl")
include("shape.jl")
include("scene.jl")
include("intersection.jl")
include("render.jl")
include("interactive.jl")

end # Rays
