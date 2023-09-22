module Rays

using Combinatorics: combinations
using IterTools: product
using KernelAbstractions
using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!
using ProgressMeter
using StatsBase: countmap

include("camera.jl")
include("shape.jl")
include("render.jl")

end # Rays
