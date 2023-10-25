module Rays

import SimpleDirectMediaLayer as SDL

using FunctionWrappers
import FunctionWrappers: FunctionWrapper

using Base.Threads
using Combinatorics: combinations
using IterTools: product
using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!
using SimpleDirectMediaLayer.LibSDL2


include("utils.jl")
include("camera.jl")
include("shape.jl")
include("scene.jl")
include("intersection.jl")
include("render.jl")
include("interactive.jl")

end # Rays
