module Rays

import SimpleDirectMediaLayer as SDL

using Accessors: @set
using Base.Threads: nthreads, threadid
using Polyester
using FunctionWrappers: FunctionWrapper
using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!, mul!
using SimpleDirectMediaLayer.LibSDL2

include("utils.jl")
include("camera.jl")
include("transform.jl")
include("shape.jl")
include("intersection.jl")
include("texture.jl")
include("scene.jl")
include("render.jl")
include("interactive.jl")

end # Rays
