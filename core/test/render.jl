using Test
using Rays: Rays
using Random: seed!
using DelimitedFiles
using LinearAlgebra: normalize, norm

@testset "blur_kernel" begin
	kernel, Δi_max = Rays.get_blur_kernel(1.6)
	@test sum(kernel) ≈ 1.0
	@test all(kernel .>= 0.0)
	@test kernel ≈ Float32[0.00029104948, 0.2252549, 0.5489081, 0.2252549, 0.00029104948]
end

@testset "shape_renders" begin
	scene = Rays.Scene()

	camera = Rays.Camera(; screen_res = [25, 25])
	from = Float32[2.0, 2.0, 2.0]
	to = zeros(Float32, 3)
	Rays.look_at!(camera, from, to)
	push!(scene, camera)

	R = 1.0f0

	function my_field(loc::Vector{F})::F where {F}
		x, y, z = loc
		out = -4.0f0
		for i ∈ 0:3
			θ = convert(F, π * (i / 2 + 0.33))
			out += one(F) / sqrt((x - cos(θ))^2 + (y - sin(θ))^2 + z^2)
		end
		return out
	end

	shapes = [
		Rays.Cube(R),
		Rays.menger_sponge(R, 4),
		Rays.Sphere(R),
		Rays.Tetrahedron(R),
		Rays.sierpinski_pyramid(R, 4),
		Rays.ImplicitSurface(
			my_field,
			R_bound = 1.5f0,
			n_divisions = 50,
			tol = 1.0f-5,
			itermax = 10,
			name = :equipotential_surface,
		),
		Rays.RevolutionSurface(
			z -> z^2 + 0.1f0,
			3.0f0,
			-1.0f0,
			0.75f0;
			n_divisions = 50,
		),
	]

	push!(scene, shapes[1])
	Rays.set_dropoff_curve_default!(scene, camera)

	function simple_view!(shape)
		Rays.clear_shapes!(scene)
		push!(scene, shape)
		Rays.render!(scene)
		return nothing
	end

	# writing files: writedlm("cube_render.csv", camera.canvas[1, :, :], ',') 

	simple_view!(shapes[1])
	@test camera.canvas[1, :, :] ≈
		  readdlm(normpath(@__DIR__, "files/cube_render.csv"), ',', Float32, '\n')

	simple_view!(shapes[2])
	@test camera.canvas[1, :, :] ≈
		  readdlm(normpath(@__DIR__, "files/Menger_sponge_render.csv"), ',', Float32, '\n')

	simple_view!(shapes[3])
	@test camera.canvas[1, :, :] ≈
		  readdlm(normpath(@__DIR__, "files/sphere_render.csv"), ',', Float32, '\n')

	simple_view!(shapes[4])
	@test camera.canvas[1, :, :] ≈
		  readdlm(normpath(@__DIR__, "files/tetrahedron_render.csv"), ',', Float32, '\n')

	simple_view!(shapes[5])
	@test camera.canvas[1, :, :] ≈ readdlm(
		normpath(@__DIR__, "files/Sierpinski_pyramid_render.csv"),
		',',
		Float32,
		'\n',
	)

	simple_view!(shapes[6])
	@test camera.canvas[1, :, :] ≈ readdlm(
		normpath(@__DIR__, "files/implicit_surface_render.csv"),
		',',
		Float32,
		'\n',
	)

	simple_view!(shapes[7])
	@test camera.canvas[1, :, :] ≈ readdlm(
		normpath(@__DIR__, "files/revolution_surface_render.csv"),
		',',
		Float32,
		'\n',
	)

	# Many shapes
	seed!(314156)
	Rays.clear_shapes!(scene)
	n_cubes = 250

	for i ∈ 1:n_cubes
		center = Float32.(1.2 * (rand(3) * 2 .- 1))
		R = 0.25f0 / ((2 * norm(center))^2 + 1)
		transform = Rays.translation(center) * Rays.rotation(normalize(rand(Float32, 3)), Float32(2π) * rand(Float32))
		cube = Rays.Cube(R)
		push!(scene, cube; transform)
	end

	Rays.render!(scene)
	@test camera.canvas[1, :, :] ≈
		  readdlm(normpath(@__DIR__, "files/many_cubes_render.csv"), ',', Float32, '\n')
end

@testset "texture renders" begin
	scene = Rays.Scene()

	camera = Rays.Camera(; screen_res = [25, 25])
	from = Float32[2.0, 2.0, 2.0]
	to = zeros(Float32, 3)
	Rays.look_at!(camera, from, to)
	push!(scene, camera)

	R = 1.0f0
	push!(scene, Rays.Cube(R))

	Rays.set_dropoff_curve_default!(scene, camera)

	julia_green = Float32[0.22, 0.596, 0.149]
	julia_purple = Float32[0.584, 0.345, 0.698]
	julia_red = Float32[0.796, 0.235, 0.2]
	julia_colors = hcat(julia_green, julia_purple, julia_red)

	function colorfield!(color::Vector{Float32}, loc_int::Vector{Float32})::Nothing
		for i ∈ 1:3
			color[i] = min(2 * sqrt(2) * abs(loc_int[i]), 1.0f0)
		end
		return nothing
	end

	Rays.set_texture!(scene, :cube, Rays.IntegerMappingTexture(julia_colors, :dim))
	Rays.render!(scene)
	@test camera.canvas ≈ reshape(
		readdlm(
			normpath(@__DIR__, "files/integer_mapping_texture_render.csv"),
			',',
			Float32,
			'\n',
		),
		size(camera.canvas),
	)

	Rays.set_texture!(scene, :cube, Rays.UniformTexture(Float32[0.0, 1.0, 0.0]))
	Rays.render!(scene)
	@test camera.canvas ≈ reshape(
		readdlm(normpath(@__DIR__, "files/uniform_texture_render.csv"), ',', Float32, '\n'),
		size(camera.canvas),
	)

	Rays.set_texture!(scene, :cube, Rays.ColorFieldTexture(colorfield!))
	Rays.render!(scene)
	@test camera.canvas ≈ reshape(
		readdlm(
			normpath(@__DIR__, "files/color_field_texture_render.csv"),
			',',
			Float32,
			'\n',
		),
		size(camera.canvas),
	)
end
