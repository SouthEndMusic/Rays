import Rays
using LinearAlgebra: cross, norm

@testset "Camera" begin
    loc = zeros(3)
    dir = [0.0, 1.0, 0.0]
    up = [0.0, 0.0, 1.0]
    screen_size = [0.1, 0.1]
    screen_dist = [0.1]
    screen_res = [100, 100]
    camera = Rays.Camera(loc, dir, up, screen_size, screen_dist, screen_res)

    @test camera isa Rays.Camera
    @test camera.right == cross(camera.dir, camera.up)
    @test camera.warp! == identity
    @test string(camera) == "<Camera 'camera'>"
end

@testset "look at" begin
    camera = Rays.Camera()
    from = Float32[0.2947051, 0.45194465, 0.05835851]
    to = Float32[0.21519046, 0.39008522, 0.3857729]
    Rays.look_at!(camera, from, to)

    @test camera.loc ≈ from
    @test norm(camera.dir) ≈ 1.0
    @test camera.dir ≈ Float32[-0.2321169, -0.1805783, 0.95577884]
    @test norm(camera.up) ≈ 1.0
    @test camera.up ≈ Float32[0.7543785, 0.5868784, 0.29408634]
    @test norm(camera.right) ≈ 1.0
    @test camera.right ≈ Float32[-0.6140316, 0.7892814, 0.0]

    from[1:2] = to[1:2]

    @test_throws "In 'look_at!', the camera cannot point straight up or down." Rays.look_at!(
        camera,
        from,
        to,
    )
end

@testset "set ray" begin
    camera = Rays.Camera()
    camera.screen_res .= [100, 100]

    from = Float32[0.61499727, 0.9667763, 0.1495682]
    to = Float32[0.87458175, 0.4226909, 0.0016991668]
    Rays.look_at!(camera, from, to)

    ray::Rays.Ray{Float32} = Rays.Ray()
    Rays.set_ray!(ray, camera, (24, 68))
    @test ray isa Rays.Ray
    @test ray.loc ≈ Float32[0.6436098, 0.8657537, 0.15174258]
    @test norm(ray.dir) ≈ 1.0
    @test ray.dir ≈ Float32[0.2724516, -0.96194667, 0.020704657]
end
