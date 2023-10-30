import Rays
using LinearAlgebra: normalize!

@testset "Sphere" begin
    center = ones(3)
    R = 0.5
    sphere = Rays.Sphere(center, R)
    @test sphere isa Rays.Sphere
    @test string(sphere) == "<Sphere 'sphere'>\n"

    intersection = Rays.Intersection()
    (; ray) = intersection
    ray.loc .= 0.0
    ray.dir .= 1.0
    normalize!(ray.dir)
    closer_intersection_found = Rays.intersect!(intersection, sphere)
    @test closer_intersection_found
    @test intersection.t[1] ≈ sqrt(3) - 0.5

    Rays.reset_intersection!(intersection)
    ray.dir .= [0.0, 0.0, 1.0]
    closer_intersection_found = Rays.intersect!(intersection, sphere)
    @test !closer_intersection_found
    @test intersection.t[1] == Inf
end

@testset "Cube" begin
    center = [0.0, 2.0, 0.0]
    R = 0.5
    cube = Rays.Cube(center, R)
    @test cube isa Rays.Cube

    cube.name[1] = :my_awesome_cube
    @test string(cube) == "<Cube 'my_awesome_cube'>\n"

    intersection = Rays.Intersection()
    (; ray) = intersection
    ray.loc .= 0.0
    ray.dir .= [0.0, 1.0, 0.0]
    normalize!(ray.dir)

    closer_intersection_found = Rays.intersect!(intersection, cube)
    @test closer_intersection_found
    @test intersection.t[1] ≈ 1.5
    @test intersection.dim[1] == 2
end
