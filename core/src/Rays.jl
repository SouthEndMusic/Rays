module Rays

import SimpleDirectMediaLayer as SDL

using FunctionWrappers
import FunctionWrappers: FunctionWrapper

using Accessors: @set
using Base.Threads
using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!
using SimpleDirectMediaLayer.LibSDL2

const ScalarFunc = FunctionWrapper{F,Tuple{F}} where {F<:AbstractFloat}
const Transform = FunctionWrapper{Nothing,Tuple{Vector{F}}} where {F<:AbstractFloat}
const ScalarField = FunctionWrapper{F,Tuple{Vector{F}}} where {F}
const VectorField = FunctionWrapper{Nothing,Tuple{Vector{F},Vector{F}}} where {F}

include("utils.jl")
include("camera.jl")
include("shape.jl")
include("Texture.jl")
include("scene.jl")
include("intersection.jl")
include("render.jl")
include("interactive.jl")

end # Rays
