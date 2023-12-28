using Test: @testset
using SafeTestsets: @safetestset

@testset "Rays" begin
    @safetestset "camera" include("camera_test.jl")
    @safetestset "shape" include("shape_test.jl")
    @safetestset "scene" include("scene_test.jl")
    @safetestset "render" include("render_test.jl")
end
