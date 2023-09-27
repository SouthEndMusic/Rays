using Test
using Rays
using LinearAlgebra: normalize!

@testset "Sphere" begin
    center = ones(3)
    R = 0.5
    sphere = Rays.Sphere(center, R)
    @test sphere isa Rays.Sphere

    loc = zeros(3)
    dir = ones(3)
    normalize!(dir)
    ray = Rays.Ray(loc, dir)
    closer_intersection_found, intersection = Rays.intersect!(ray, sphere)
    @test closer_intersection_found
    @test intersection.t[1] ≈ sqrt(3) - 0.5

    ray.dir[:] = [0.0, 0.0, 1.0]
    closer_intersection_found, intersection = Rays.intersect!(ray, sphere)
    @test !closer_intersection_found
    @test intersection.t[1] == Inf
end

@testset "Cube" begin
    center = [0.0, 2.0, 0.0]
    R = 0.5
    cube = Rays.Cube(center, R)
    @test cube isa Rays.Cube

    loc = zeros(3)
    dir = [0.0, 1.0, 0.0]
    normalize!(dir)
    ray = Rays.Ray(loc, dir)

    closer_intersection_found, intersection = Rays.intersect!(ray, cube)
    @test closer_intersection_found
    @test intersection.t[1] ≈ 1.5
    @test intersection.dim[1] == 2
end
