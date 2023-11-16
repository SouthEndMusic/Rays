using Test
using Rays: Rays
using Accessors: @set
using LinearAlgebra: normalize!, norm

@testset "Sphere" begin
	center = ones(3)
	R = 0.5
	sphere = Rays.Sphere(center, R)
	@test sphere isa Rays.Sphere
	@test string(sphere) == "<Sphere 'sphere'>"

	intersection = Rays.Intersection()
	(; ray) = intersection
	ray.loc .= 0.0
	ray.dir .= 1.0
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

@testset "Cube" begin
	center = [0.0, 2.0, 0.0]
	R = 0.5
	cube = Rays.Cube(center, R)
	cube = @set cube.name = :my_awesome_cube

	@test cube isa Rays.Cube
	@test string(cube) == "<Cube 'my_awesome_cube'>"

	intersection = Rays.Intersection()
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
	center = zeros(3)
	R = 0.5
	sponge = Rays.menger_sponge(center, R, 4)
	sponge = @set sponge.name = :my_awesome_sponge

	@test sponge isa Rays.FractalShape
	@test string(sponge) == "<FractalShape 'my_awesome_sponge'; 20 subshapes of type Rays.Cube{Float64}>"

	intersection = Rays.Intersection()
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
	center = zeros(3)
	R = 0.5
	function f(x)
		return norm(x) - R
	end
	sphere = Rays.ImplicitSurface(f, center, name = :sphere, R_bound = 1.1 * R, tol = 1e-5)
	@test sphere isa Rays.ImplicitSurface
	@test string(sphere) == "<ImplicitSurface 'sphere'; function 'f' and finite difference gradient>"

	intersection = Rays.Intersection()
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
	center = zeros(3)
	vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
	faces = [1 2 3;]
	triangle = Rays.TriangleShape(vertices, faces, center; name = :my_awesome_triangle)
	@test string(triangle) == "<TriangleShape 'my_awesome_triangle'; with 3 vertices and 1 faces>"

	intersection = Rays.Intersection()
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
