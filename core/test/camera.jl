using Test
using Rays
using LinearAlgebra: cross, norm

@testset "Camera" begin
    loc = zeros(3)
    dir = [0.0, 1.0, 0.0]
    up = [0.0, 0.0, 1.0]
    screen_size = [0.1, 0.1]
    screen_dist = [0.1]
    screen_res = [100, 100]
    camera = Rays.Camera([loc], [dir], [up], [screen_size], [screen_dist], [screen_res])

    @test camera isa Rays.Camera
    @test all(@. camera.right == cross(camera.dir, camera.up))
    @test camera.warp == fill(identity, length(camera.loc))
end

@testset "look at" begin
    camera = Rays.Camera(1)
    from = [0.2947051, 0.45194465, 0.05835851]
    to = [0.21519046, 0.39008522, 0.3857729]
    Rays.look_at!(camera, 1, from, to)

    @test camera.loc[1] == from
    @test norm(camera.dir[1]) ≈ 1.0
    @test camera.dir[1] ≈ Float32[-0.2321169, -0.1805783, 0.95577884]
    @test norm(camera.up[1]) ≈ 1.0
    @test camera.up[1] ≈ Float32[0.7543785, 0.5868784, 0.29408634]
    @test norm(camera.right[1]) ≈ 1.0
    @test camera.right[1] ≈ Float32[-0.6140316, 0.7892814, 0.0]

    from[1:2] = to[1:2]

    @test_throws "In 'look_at!', the camera cannot point straight up or down." Rays.look_at!(
        camera,
        1,
        from,
        to,
    )
end

@testset "get ray" begin
    camera = Rays.Camera(1)
    camera.screen_res[1] = [100, 100]

    from = [0.61499727, 0.9667763, 0.1495682]
    to = [0.87458175, 0.4226909, 0.0016991668]
    Rays.look_at!(camera, 1, from, to)

    ray = Rays.get_ray(camera, 1, [24, 68])
    @test ray isa Rays.Ray
    @test ray.loc ≈ Float32[0.6827902, 0.8868463, 0.1429134]
    @test norm(ray.dir) ≈ 1.0
    @test ray.dir ≈ Float32[0.64553064, -0.761101, -0.0633676]

    @test_throws "Pixel indices must fall inside screen res [100, 100], got [0, 68]." Rays.get_ray(
        camera,
        1,
        [0, 68],
    )

    @test_throws "Pixel indices must fall inside screen res [100, 100], got [123, 68]." Rays.get_ray(
        camera,
        1,
        [123, 68],
    )
end
