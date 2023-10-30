using Test: @testset
using SafeTestsets: @safetestset

@testset "Rays" begin
    @safetestset "camera" include("camera.jl")
    @safetestset "shape" include("shape.jl")
    @safetestset "scene" include("scene.jl")
    @safetestset "render" include("render.jl")
end
