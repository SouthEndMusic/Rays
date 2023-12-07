using Test
using Rays: Rays
using Accessors: @set
using LinearAlgebra: normalize!, norm

@testset "Sphere" begin
	R = 0.5
	sphere = Rays.Sphere(R)
	@test sphere isa Rays.Sphere
	@test string(sphere) == "<Sphere 'sphere'>"

	matrix_prototype = zeros(3, 3)
	intersection = Rays.Intersection(1; matrix_prototype)
	ray_loc,
	ray_dir,
	ray_camera_loc,
	ray_camera_dir,
	t,
	cache_int,
	cache_float,
	color,
	name_intersected = Rays.get_caches(intersection, 1)


	ray_camera_loc .= 0.0
	ray_camera_dir .= 1.0
	normalize!(ray_camera_dir)
	translation = Rays.translation(ones(3))
	intersector! = Rays.create_intersector(sphere, translation; matrix_prototype)
	intersector!(
		t,
		ray_loc,
		ray_dir,
		ray_camera_loc,
		ray_camera_dir,
		cache_int,
		cache_float,
		name_intersected,
	)
	@test name_intersected[1] == :sphere
	@test t[1] ≈ sqrt(3) - 0.5

	Rays.reset_intersection!(t, name_intersected)
	ray_camera_dir .= [0.0, 0.0, 1.0]
	intersector!(t,
		ray_loc,
		ray_dir,
		ray_camera_loc,
		ray_camera_dir,
		cache_int,
		cache_float,
		name_intersected,
	)
	@test name_intersected[1] == :none
	@test t[1] == Inf
end

@testset "Cube" begin
	R = 0.5
	cube = Rays.Cube(R)
	cube = @set cube.name = :my_awesome_cube

	@test cube isa Rays.Cube
	@test string(cube) == "<Cube 'my_awesome_cube'>"

	matrix_prototype = zeros(3, 3)
	intersection = Rays.Intersection(1; matrix_prototype)
	ray_loc,
	ray_dir,
	ray_camera_loc,
	ray_camera_dir,
	t,
	cache_int,
	cache_float,
	color,
	name_intersected = Rays.get_caches(intersection, 1)

	ray_camera_loc .= 0.0
	ray_camera_dir .= [0.0, 1.0, 0.0]
	normalize!(ray_camera_dir)
	translation = Rays.translation([0.0, 2.0, 0.0])
	intersector! = Rays.create_intersector(cube, translation; matrix_prototype)

	intersector!(
		t,
		ray_loc,
		ray_dir,
		ray_camera_loc,
		ray_camera_dir,
		cache_int,
		cache_float,
		name_intersected,
	)

	@test name_intersected[1] == :my_awesome_cube
	@test t[1] ≈ 1.5
	@test cache_int[1] == 2
end

@testset "FractalShape" begin
	R = 0.5
	sponge = Rays.menger_sponge(R, 4)
	sponge = @set sponge.name = :my_awesome_sponge

	@test sponge isa Rays.FractalShape
	@test string(sponge) ==
		  "<FractalShape 'my_awesome_sponge'; 20 subshapes of <Cube 'cube'>>"

	matrix_prototype = zeros(3, 3)
	intersection = Rays.Intersection(1; matrix_prototype)
	ray_loc,
	ray_dir,
	ray_camera_loc,
	ray_camera_dir,
	t,
	cache_int,
	cache_float,
	color,
	name_intersected = Rays.get_caches(intersection, 1)

	ray_loc .= [1.0, 2.0, 3.0]
	ray_dir .= -ray_loc
	normalize!(ray_dir)

	closer_intersection_found = Rays._intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, sponge)
	@test closer_intersection_found
	@test !isinf(t[1])
	@test cache_int[1] == 3
end

@testset "ImplicitSurface" begin
	R = 0.5
	function f(x)
		return norm(x) - R
	end
	vector_prototype = zeros(3)
	sphere = Rays.ImplicitSurface(f, name = :sphere, R_bound = 1.1 * R, tol = 1e-5; vector_prototype)
	@test sphere isa Rays.ImplicitSurface
	@test string(sphere) ==
		  "<ImplicitSurface 'sphere'; function 'f' and finite difference gradient>"

	matrix_prototype = zeros(3, 3)
	intersection = Rays.Intersection(1; matrix_prototype)
	ray_loc,
	ray_dir,
	ray_camera_loc,
	ray_camera_dir,
	t,
	cache_int,
	cache_float,
	color,
	name_intersected = Rays.get_caches(intersection, 1)

	ray_loc .= 1.0
	ray_dir .= -1.0
	normalize!(ray_dir)
	closer_intersection_found = Rays._intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, sphere)
	@test closer_intersection_found
	@test t[1] ≈ sqrt(3) - 0.5

	Rays.reset_intersection!(t, name_intersected)
	ray_dir .= [0.0, 0.0, 1.0]
	closer_intersection_found = Rays._intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, sphere)
	@test !closer_intersection_found
	@test t[1] == Inf
end

@testset "TriangleShape" begin
	vertices = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
	faces = [1 2 3;]
	triangle = Rays.TriangleShape(vertices, faces; name = :my_awesome_triangle)
	@test string(triangle) ==
		  "<TriangleShape 'my_awesome_triangle'; with 3 vertices and 1 faces>"

	matrix_prototype = zeros(3, 3)
	intersection = Rays.Intersection(1; matrix_prototype)
	ray_loc,
	ray_dir,
	ray_camera_loc,
	ray_camera_dir,
	t,
	cache_int,
	cache_float,
	color,
	name_intersected = Rays.get_caches(intersection, 1)

	ray_loc .= [0.25, 0.25, 1.0]
	ray_dir .= [0.0, 0.0, -1.0]
	closer_intersection_found = Rays._intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, triangle)
	@test closer_intersection_found
	@test t[1] ≈ 1.0

	Rays.reset_intersection!(t, name_intersected)
	ray_loc .= [0.75, 0.75, 1.0]
	closer_intersection_found = Rays._intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, triangle)
	@test !closer_intersection_found
end

@testset "RevolutionSurface" begin
	r(z) = z^2
	revolution_shape =
		Rays.RevolutionSurface(r, 1.0, -1.0, 1.0; tol = 1e-7, n_divisions = 50)
	@test string(revolution_shape) ==
		  "<RevolutionSurface 'revolution_surface'; function 'r' and finite difference derivative>"

	matrix_prototype = zeros(3, 3)
	intersection = Rays.Intersection(1; matrix_prototype)
	ray_loc,
	ray_dir,
	ray_camera_loc,
	ray_camera_dir,
	t,
	cache_int,
	cache_float,
	color,
	name_intersected = Rays.get_caches(intersection, 1)

	ray_loc .= [1.0, 1.0, 0.5]
	ray_dir .= [-1.0, -1.0, 1e-8]
	normalize!(ray_dir)
	closer_intersection_found = Rays._intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, revolution_shape)
	@test closer_intersection_found
	@test t[1] ≈ sqrt(2.0) - 0.25
end
