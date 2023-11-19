using Test
using Rays: Rays
using Accessors: @set
using LinearAlgebra: normalize!, norm

@testset "Sphere" begin
    R = 0.5
    sphere = Rays.Sphere(R)
    @test sphere isa Rays.Sphere
    @test string(sphere) == "<Sphere 'sphere'>"

    intersection = Rays.Intersection(; F = Float64)
    (; ray_camera) = intersection
    ray_camera.loc .= 0.0
    ray_camera.dir .= 1.0
    normalize!(ray_camera.dir)
    translation = Rays.get_translation(ones(3))
    intersector! = Rays.create_intersector(sphere, translation)
    intersector!(intersection)
    @test intersection.name_intersected[1] == :sphere
    @test intersection.t[1] ≈ sqrt(3) - 0.5

    Rays.reset_intersection!(intersection)
    ray_camera.dir .= [0.0, 0.0, 1.0]
    intersector!(intersection)
    @test intersection.name_intersected[1] == :none
    @test intersection.t[1] == Inf
end

@testset "Cube" begin
    R = 0.5
    cube = Rays.Cube(R)
    cube = @set cube.name = :my_awesome_cube

    @test cube isa Rays.Cube
    @test string(cube) == "<Cube 'my_awesome_cube'>"

    intersection = Rays.Intersection(; F = Float64)
    (; ray) = intersection
    ray.loc .= 0.0
    ray.dir .= [0.0, 1.0, 0.0]
    normalize!(ray.dir)

    closer_intersection_found = Rays._intersect_ray!(intersection, cube)
    @test closer_intersection_found
    @test intersection.t[1] ≈ 1.5
    @test intersection.dim[1] == 2
end

@testset "FractalShape" begin
    R = 0.5
    sponge = Rays.menger_sponge(R, 4)
    sponge = @set sponge.name = :my_awesome_sponge

    @test sponge isa Rays.FractalShape
    @test string(sponge) ==
          "<FractalShape 'my_awesome_sponge'; 20 subshapes of type Rays.Cube{Float64}>"

    intersection = Rays.Intersection(; F = Float64)
    (; ray) = intersection
    ray.loc .= [1.0, 2.0, 3.0]
    ray.dir .= -ray.loc
    normalize!(ray.dir)

    closer_intersection_found = Rays._intersect_ray!(intersection, sponge)
    @test closer_intersection_found
    @test !isinf(intersection.t[1])
    @test intersection.dim[1] == 2
end

@testset "ImplicitSurface" begin
    R = 0.5
    function f(x)
        return norm(x) - R
    end
    sphere = Rays.ImplicitSurface(f, name = :sphere, R_bound = 1.1 * R, tol = 1e-5)
    @test sphere isa Rays.ImplicitSurface
    @test string(sphere) ==
          "<ImplicitSurface 'sphere'; function 'f' and finite difference gradient>"

    intersection = Rays.Intersection(; F = Float64)
    (; ray) = intersection
    ray.loc .= 1.0
    ray.dir .= -1.0
    normalize!(ray.dir)
    closer_intersection_found = Rays._intersect_ray!(intersection, sphere)
    @test closer_intersection_found
    @test intersection.t[1] ≈ sqrt(3) - 0.5

    Rays.reset_intersection!(intersection)
    ray.dir .= [0.0, 0.0, 1.0]
    closer_intersection_found = Rays._intersect_ray!(intersection, sphere)
    @test !closer_intersection_found
    @test intersection.t[1] == Inf
end

@testset "TriangleShape" begin
    vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
    faces = [1 2 3;]
    triangle = Rays.TriangleShape(vertices, faces; name = :my_awesome_triangle)
    @test string(triangle) ==
          "<TriangleShape 'my_awesome_triangle'; with 3 vertices and 1 faces>"

    intersection = Rays.Intersection(; F = Float64)
    (; ray) = intersection
    ray.loc .= [0.25, 0.25, 1.0]
    ray.dir .= [0.0, 0.0, -1.0]
    closer_intersection_found = Rays._intersect_ray!(intersection, triangle)
    @test closer_intersection_found
    @test intersection.t[1] ≈ 1.0

    Rays.reset_intersection!(intersection)
    ray.loc .= [0.75, 0.75, 1.0]
    closer_intersection_found = Rays._intersect_ray!(intersection, triangle)
    @test !closer_intersection_found
end

@testset "RevolutionSurface" begin
    r(z) = z^2
    revolution_shape =
        Rays.RevolutionSurface(r, 1.0, -1.0, 1.0; tol = 1e-7, n_divisions = 50)
    @test string(revolution_shape) ==
          "<RevolutionSurface 'revolution_surface'; function 'r' and finite difference derivative>"

    intersection = Rays.Intersection(; F = Float64)
    (; ray) = intersection
    ray.loc .= [1.0, 1.0, 0.5]
    ray.dir .= [-1.0, -1.0, 1e-8]
    normalize!(ray.dir)
    closer_intersection_found = Rays._intersect_ray!(intersection, revolution_shape)
    @test closer_intersection_found
    @test intersection.t[1] ≈ sqrt(2.0) - 0.25
end
