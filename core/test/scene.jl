using Test
using Rays: Rays

@testset "Scene string representation" begin
    scene = Rays.Scene()
    push!(scene, Rays.Camera())
    push!(scene, Rays.Cube(1.0f0))
    Rays.partition!(scene)
    @test string(scene) ==
          "Scene (array type Matrix{Float32}):\n* Partition:\n\t<Partition;\n\t total number of objects: 1,\n\t number of leaf nodes: 2,\n\t average number of objects per leaf node: 1,\n\t maximum number of objects per leaf node: 1,\n\t number of leaf nodes per depth:\n\t\t1: 2\n>\n\n* Cameras:\n\t<Camera 'camera'>\n\n* Shapes:\n\t<Cube 'cube'>\n\n* Transforms:\n\tcube: <AffineTransform; identity>\n\n* Textures:\n\tcube: <UniformTexture; color = Float32[1.0, 1.0, 1.0]>\n"
end
