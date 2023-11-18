using Test
using Rays: Rays
using Random: seed!
using DelimitedFiles

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

    origin = zeros(Float32, 3)
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
        Rays.Cube(origin, R),
        Rays.menger_sponge(origin, R, 4),
        Rays.Sphere(origin, R),
        Rays.Tetrahedron(origin, R),
        Rays.sierpinski_pyramid(origin, R, 4),
        Rays.ImplicitSurface(
            my_field,
            origin,
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
            0.75f0,
            zeros(Float32, 3);
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

    simple_view!(shapes[1])
    @test camera.canvas[1, :, :] ≈ readdlm("files/cube_render.csv", ',', Float32, '\n')

    simple_view!(shapes[2])
    @test camera.canvas[1, :, :] ≈
          readdlm("files/Menger_sponge_render.csv", ',', Float32, '\n')

    simple_view!(shapes[3])
    @test camera.canvas[1, :, :] ≈ readdlm("files/sphere_render.csv", ',', Float32, '\n')

    simple_view!(shapes[4])
    @test camera.canvas[1, :, :] ≈
          readdlm("files/tetrahedron_render.csv", ',', Float32, '\n')

    simple_view!(shapes[5])
    @test camera.canvas[1, :, :] ≈
          readdlm("files/Sierpinski_pyramid_render.csv", ',', Float32, '\n')

    simple_view!(shapes[6])
    @test camera.canvas[1, :, :] ≈
          readdlm("files/implicit_surface_render.csv", ',', Float32, '\n')

    simple_view!(shapes[7])
    @test camera.canvas[1, :, :] ≈
          readdlm("files/revolution_surface_render.csv", ',', Float32, '\n')

    # Many shapes
    seed!(31415)
    Rays.clear_shapes!(scene)
    n_cubes = 250

    for i ∈ 1:n_cubes
        center = rand(Float32, 3) * 2 .- 1
        R = rand(Float32) / 10
        cube = Rays.Cube(center, R)
        push!(scene, cube)
    end

    Rays.render!(scene)
    @test camera.canvas[1, :, :] ≈
          readdlm("files/many_cubes_render.csv", ',', Float32, '\n')
end
