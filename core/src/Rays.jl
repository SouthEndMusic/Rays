module Rays

using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!
using Images
using VideoIO
using ProgressMeter
using StatsBase: countmap

include("camera.jl")
include("shape.jl")
include("render.jl")

end # Rays
