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

end

@testset "Menger sponge" begin

end
