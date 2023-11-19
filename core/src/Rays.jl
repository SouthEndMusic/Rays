module Rays

import SimpleDirectMediaLayer as SDL

using Accessors: @set
using Base.Threads
using FunctionWrappers: FunctionWrapper
using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!
using SimpleDirectMediaLayer.LibSDL2

const ScalarFunc = FunctionWrapper{F,Tuple{F}} where {F<:AbstractFloat}
const Transform = FunctionWrapper{Nothing,Tuple{Vector{F}}} where {F<:AbstractFloat}
const ScalarField = FunctionWrapper{F,Tuple{Vector{F}}} where {F<:AbstractFloat}
const VectorField =
    FunctionWrapper{Nothing,Tuple{Vector{F},Vector{F}}} where {F<:AbstractFloat}

include("utils.jl")
include("camera.jl")
include("transforms.jl")
include("shape.jl")
include("intersection.jl")
include("texture.jl")
include("scene.jl")
include("render.jl")
include("interactive.jl")

end # Rays
