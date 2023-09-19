using Rays
using Images
using ProgressMeter
using VideoIO

function Sierpinski_pyramid_view()
    julia_blue = [0.251, 0.388, 0.847]
    julia_green = [0.22, 0.596, 0.149]
    julia_purple = [0.584, 0.345, 0.698]
    julia_red = [0.796, 0.235, 0.2]

    julia_colors = [julia_blue, julia_green, julia_purple, julia_red]

    camera = Rays.Camera(1)
    camera.screen_res[1] = [1000, 1000]
    camera.screen_dist[1] = [1.0]
    from = [-5.0, -5.0, 5.0]
    to = zeros(3)
    Rays.look_at!(camera, 1, from, to)
    dropoff_curve(t) = clamp(7.0 - 0.9 * (t - 1), 0, 1.0)

    tetrahedron = Rays.Sierpinski_pyramid(zeros(3), 0.5, 4)
    collect_metadata = Dict(:face_int => Int)
    t_int, metadata = Rays.shape_view(camera, 1, tetrahedron; collect_metadata)

    screen_res = camera.screen_res[1]
    color = zeros(3, screen_res...)
    face_int = metadata[:face_int]
    for i = 1:screen_res[1]
        for j = 1:screen_res[2]
            face_int_ = face_int[i, j]
            color[:, i, j] = if iszero(face_int_)
                zeros(3)
            else
                julia_colors[face_int_]
            end
        end
    end

    canvas = Rays.cam_is_source(t_int, dropoff_curve; color)
    canvas = permutedims(canvas, [1, 3, 2])
    reverse!(canvas)
    colorview(RGB, canvas)
end

# Simple cube and sphere view
function simple_view()
    julia_green = [0.22, 0.596, 0.149]
    julia_purple = [0.584, 0.345, 0.698]
    julia_red = [0.796, 0.235, 0.2]

    julia_colors = [julia_green, julia_purple, julia_red]

    loc = zeros(3)
    dir = zeros(3)
    up = zeros(3)
    screen_size = [0.1, 0.1]
    screen_dist = [0.2]
    screen_res = [450, 450]
    dropoff_curve(t) = clamp(1 - 0.3 * (t - 1), 0, 1)
    focus_curve(t) = 0.5 + 20 * abs(t - 3)
    function warp!(v::Vector{Float64})::Nothing
        v[3] = v[3] + 0.1 * sin(250 * v[2])
        v[1] = v[1] + 0.1 * sin(250 * v[2])
        return nothing
    end
    cam = Rays.Camera(
        [loc],
        [dir],
        [up],
        [screen_size],
        [screen_dist],
        [screen_res];
        warp = [warp!],
    )
    Rays.look_at!(cam, 1, [2.0, 2.0, 2.0], zeros(3))

    sphere = Rays.Sphere(zeros(3), 1.0)
    cube = Rays.Cube(zeros(3), 0.5)
    sponge = Rays.Menger_sponge(zeros(3), 0.5, 4)

    collect_metadata = Dict(:dim_int => Int)
    t_int, metadata = Rays.shape_view(cam, 1, sponge; collect_metadata)
    dim_int = metadata[:dim_int]
    color = zeros(3, screen_res...)
    for i = 1:screen_res[1]
        for j = 1:screen_res[2]
            dim_int_ = dim_int[i, j]
            color[:, i, j] = if iszero(dim_int_)
                zeros(3)
            else
                julia_colors[dim_int_]
            end
        end
    end

    canvas = Rays.cam_is_source(t_int, dropoff_curve; color)
    canvas = Rays.add_depth_of_field(canvas, t_int, focus_curve)
    canvas = permutedims(canvas, [1, 3, 2])
    reverse!(canvas)

    colorview(RGB, canvas)
end

# Making animation
function rotation_animation()
    julia_green = [0.22, 0.596, 0.149]
    julia_purple = [0.584, 0.345, 0.698]
    julia_red = [0.796, 0.235, 0.2]

    julia_colors = [julia_green, julia_purple, julia_red]

    loc = zeros(3)
    dir = zeros(3)
    up = zeros(3)
    screen_size = [0.1, 0.1]
    screen_dist = [0.2]
    screen_res = [250, 250]
    dropoff_curve(t) = clamp(1 - 0.3 * (t - 1), 0, 1)
    focus_curve(t) = 0.5 + 20 * abs(t - 3)
    function warp!(v::Vector{Float64})::Nothing
        v[3] = v[3] + 0.1 * sin(250 * v[2])
        v[1] = v[1] + 0.1 * sin(250 * v[2])
        return nothing
    end
    cam = Rays.Camera(
        [loc],
        [dir],
        [up],
        [screen_size],
        [screen_dist],
        [screen_res];
        warp = [warp!],
    )
    sponge = Rays.Menger_sponge(zeros(3), 0.5, 4)

    n_frames = 50
    framerate = 10
    ϕ = π / 3
    R = 4

    encoder_options = (; crf = 0)

    open_video_out(
        normpath(@__DIR__, "Menger_sponge_spin.mp4"),
        RGB{N0f8},
        (screen_res...,),
        framerate = framerate,
        encoder_options = encoder_options,
    ) do writer
        @showprogress "Computing frames..." for θ in range(π / 6, 7π / 6, n_frames)
            x = cos(θ) * sin(ϕ)
            y = sin(θ) * sin(ϕ)
            z = cos(ϕ)

            from = R * [x, y, z]
            to = zeros(3)

            Rays.look_at!(cam, 1, from, to)

            collect_metadata = Dict(:dim_int => Int)
            t_int, metadata = Rays.shape_view(cam, 1, sponge; collect_metadata)
            color = zeros(3, screen_res...)
            dim_int = metadata[:dim_int]
            for i = 1:screen_res[1]
                for j = 1:screen_res[2]
                    dim_int_ = dim_int[i, j]
                    color[:, i, j] = if iszero(dim_int_)
                        zeros(3)
                    else
                        julia_colors[dim_int_]
                    end
                end
            end

            canvas = Rays.cam_is_source(t_int, dropoff_curve; color)
            canvas = permutedims(canvas, [1, 3, 2])
            reverse!(canvas)
            canvas = RGB{N0f8}.([canvas[channel, :, :] for channel = 1:3]...)
            write(writer, colorview(RGB, canvas))
        end
    end
end
