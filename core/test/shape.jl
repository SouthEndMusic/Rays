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
    t_int, _ = Rays.intersect(ray, sphere)
    @test t_int ≈ sqrt(3) - 0.5

    ray.dir[:] = [0.0, 0.0, 1.0]
    t_int, _ = Rays.intersect(ray, sphere)
    @test t_int == Inf
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

    t_int, int_metadata = Rays.intersect(ray, cube)
    @test t_int ≈ 1.5
    @test haskey(int_metadata, :dim_int)
    @test int_metadata.dim_int == 2
end
