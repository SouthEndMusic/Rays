using Test
using Rays: Rays

@testset "Scene string representation" begin
    scene = Rays.Scene()
    push!(scene, Rays.Camera())
    push!(scene, Rays.Cube(zeros(Float32, 3), 1.0f0))
    @test string(scene) ==
          "Rays.Scene{Float32}\n* Cameras:\n\t<Camera 'camera'>\n\n* Shapes:\n\t<Cube 'cube'>\n"
end
