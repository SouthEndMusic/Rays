using LinearAlgebra: cross, normalize, normalize!, norm, dot, transpose!
using Images
using VideoIO
using ProgressMeter
using StatsBase: countmap

"""
Struct for possibly multiple cameras. Per camera:
loc: The location of the camera in space
dir: The direction the camera is pointing (unit length)
up: The upwards direction 
right: The right direction
screen_size: The size of the screen in world units
screen_dist: The distance between loc and the image plane
screen_res: The resolution of the resulting render
warp: A function which changes the origin of rays as they leave the camera

dir, up and right are orthormal and completely fix the location and
position of the camera.
"""
struct Camera
    loc::Vector{Vector{Float64}}
    dir::Vector{Vector{Float64}}
    up::Vector{Vector{Float64}}
    right::Vector{Vector{Float64}}
    screen_size::Vector{Vector{Float64}} # x,y
    screen_dist::Vector{Vector{Float64}}
    screen_res::Vector{Vector{Int}} # x,y
    warp::Vector{Function}
end

"""
Construct a camera object where the right vectors are computed using the cross product
and warp functions are optional.
"""
function Camera(loc, dir, up, screen_size, screen_dist, screen_res; warp = nothing)::Camera
    right = [cross(dir_, up_) for (dir_, up_) in zip(dir, up)]

    if isnothing(warp)
        warp = fill(identity, length(loc))
    end

    return Camera(loc, dir, up, right, screen_size, screen_dist, screen_res, warp)
end

"""
Point the camera with the given cam_index from 'from' to 'to'. 
'up' is defined to be a linear combination of e_z = [0,0,1] and dir
of unit length and orthogonal to dir. This does not work if dir and e_z are proportional,
e.g. the camera points straight up or straight down.
"""
function look_at!(
    camera::Camera,
    cam_index::Int,
    from::Vector{Float64},
    to::Vector{Float64},
)::Nothing
    (; loc, dir, up, right) = camera

    loc[cam_index] = from

    dir_ = normalize(to - from)
    dir_z = dir_[3]
    denom = sqrt(1 - dir_z^2)

    if denom ≈ 0.0
        error("In 'look_at!', the camera cannot point straight up or down.")
    end

    dir[cam_index] = dir_
    up[cam_index] = -dir_z * dir_
    up[cam_index][3] += 1
    up[cam_index] /= denom
    right[cam_index] = cross(dir_, up[cam_index])

    return nothing
end

"""
Get a Matrix{Float64} with the resolution of the provided camera.
"""
get_canvas(camera::Camera, cam_index::Int) = zeros(Float64, camera.screen_res[cam_index]...)

"""
loc: The origin of the ray 
dir: The direction of the ray (unit vector)
"""
struct Ray
    loc::Vector{Float64}
    dir::Vector{Float64}
end

abstract type Shape end

struct Sphere <: Shape
    center::Vector{Float64}
    R::Float64
    Rsq::Float64
end

Sphere(center::Vector{Float64}, R::Float64) = Sphere(center, R, R^2)

struct Cube <: Shape
    center::Vector{Float64}
    R::Float64
end

"""
The Menger sponge struct is similar to the Cube struct but has 2 extra fields:
depth: The recursion depth of the Menger sponge fractal
cubes: The subcubes the largest cube consists of, e.g. 1 recursion step
"""
struct Menger_sponge <: Shape
    center::Vector{Float64}
    R::Float64
    depth::Int
    cubes::Vector{Cube}
end

"""
Create a Menger sponge of given location, size and recursion depth
with the 'cubes' array automatically generated.
"""
function Menger_sponge(center::Vector{Float64}, R::Float64, depth::Int)::Menger_sponge
    cubes = Cube[]
    R_cube = R / 3

    for i = 1:3
        for j = 1:3
            for k = 1:3
                m = countmap([i, j, k])
                if 2 ∈ keys(m) && countmap([i, j, k])[2] > 1
                    continue
                end
                center_cube = zeros(3)
                center_cube += center
                center_cube += @. ([i, j, k] - 2) * 2 * R_cube

                cube = Cube(center_cube, R_cube)
                push!(cubes, cube)
            end
        end
    end
    return Menger_sponge(center, R, depth, cubes)
end

"""
Compute the location and direction of a ray emited from the camera for 
the given pixel indices.
If a warp function is provided, this function is applied to the ray location.
"""
function get_ray(camera::Camera, cam_index::Int, pixel_indices::Vector{Int})::Ray
    (; loc, dir, up, right, screen_size, screen_dist, screen_res, warp) = camera

    loc = loc[cam_index]
    dir = dir[cam_index]
    up = up[cam_index]
    right = right[cam_index]
    screen_size = screen_size[cam_index]
    screen_dist = screen_dist[cam_index][1]
    screen_res = screen_res[cam_index]
    warp! = warp[cam_index]

    if !all(pixel_indices .<= screen_res)
        error("Pixel indices must fall inside screen res $screen_res, got $pixel_indices.")
    end

    i, j = pixel_indices
    s_w, s_h = screen_size
    res_w, res_h = screen_res

    ray_loc = loc + screen_dist * dir
    ray_loc += ((j - 1) / (res_h - 1) - 0.5) * s_h * up
    ray_loc += ((i - 1) / (res_w - 1) - 0.5) * s_w * right

    ray_dir = ray_loc - loc
    normalize!(ray_dir)

    warp!(ray_loc)

    return Ray(ray_loc, ray_dir)
end

"""
Compute the intersection of a ray with a sphere as the smallest
real solution to a quadratic polynomial, if it exists.
"""
function intersect(ray::Ray, sphere::Sphere)::Tuple{Float64,NamedTuple}
    (; loc, dir) = ray
    (; center, Rsq) = sphere

    diff = loc - center
    a = norm(dir)^2
    b = 2 * dot(dir, diff)
    c = norm(diff)^2 - Rsq
    discr = b^2 - 4 * a * c

    t_int = if discr >= 0
        t_int = (-b - sqrt(discr)) / (2 * a)
    else
        Inf
    end

    return t_int, (;)
end

"""
Compute the intersection of a ray with a cube by computing the intersections
with each of the 6 face planes and then checking whether the intersection is within the face.
"""
function intersect(ray::Ray, cube::Cube)::Tuple{Float64,NamedTuple}
    (; loc, dir) = ray
    (; center, R) = cube

    diff = center - loc

    t_int = Inf
    dim_int = 0

    for r in [-R, R]
        for dim = 1:3
            t_int_candidate = (diff[dim] + r) / dir[dim]
            valid_candidate = true
            for dim_other = 1:3
                if dim == dim_other
                    continue
                end
                loc_intersect = loc[dim_other] + t_int_candidate * dir[dim_other]
                if abs(loc_intersect - center[dim_other]) > R
                    valid_candidate = false
                    break
                end
            end
            if valid_candidate
                if t_int_candidate < t_int
                    dim_int = dim
                    t_int = t_int_candidate
                end
            end
        end
    end
    return t_int, (; dim_int)
end

"""
Compute the intersection of a ray with a Menger sponge.
This is done recursively until the recursion depth of the Menger sponge.
To compute the intersection of a ray with a sub-cube, the ray location is transformed.
"""
function intersect(
    ray::Ray,
    menger_sponge::Menger_sponge;
    current_depth::Int = 0,
)::Tuple{Float64,NamedTuple}
    t_int = Inf
    dim_int = 0
    cube_intersect = nothing

    for cube in menger_sponge.cubes
        t_int_candidate, int_metadata = intersect(ray, cube)

        if t_int_candidate < t_int
            cube_intersect = cube
            t_int = t_int_candidate
            dim_int = int_metadata.dim_int
        end
    end

    if !isnothing(cube_intersect) && current_depth < menger_sponge.depth
        t_int = Inf
        loc_transformed = 3 * (ray.loc + menger_sponge.center - cube_intersect.center)
        ray_transformed = Ray(loc_transformed, ray.dir)
        t_int_candidate, int_metadata =
            intersect(ray_transformed, menger_sponge; current_depth = current_depth + 1)
        t_int_candidate /= 3.0
        if t_int_candidate < t_int
            dim_int = int_metadata.dim_int
            t_int = t_int_candidate
        end
    end

    return t_int, (; dim_int)
end

function shape_view(
    camera::Camera,
    cam_index::Int,
    shape::Shape;
    collect_metadata::Dict{Symbol,DataType} = Dict{Symbol,DataType}(),
)::Tuple{Matrix{Float64},Dict{Symbol,Matrix}}

    (; screen_res) = camera
    screen_res = screen_res[cam_index]

    metadata = Dict{Symbol,Matrix}()
    for (s, T) in collect_metadata
        metadata[s] = zeros(T, screen_res...)
    end

    t_int = get_canvas(camera, cam_index)

    for i = 1:screen_res[1]
        for j = 1:screen_res[2]
            ray = get_ray(camera, cam_index, [i, j])
            t_int_, int_metadata = intersect(ray, shape)
            t_int[i, j] = t_int_

            for s in keys(collect_metadata)
                metadata[s][i, j] = getfield(int_metadata, s)
            end
        end
    end

    return t_int, metadata
end

function cam_is_source(intersections::Matrix{Float64}, dropoff_curve; color = nothing)
    where_intersect = .!isinf.(intersections)
    t_intersect = intersections[where_intersect]
    t_min = minimum(t_intersect)
    t_max = maximum(t_intersect)


    canvas = zeros(Float64, size(intersections)...)
    canvas[where_intersect] = @. dropoff_curve(t_intersect)

    if !isnothing(color)
        canvas_color = zeros(Float64, 3, size(intersections)...)

        for channel = 1:3
            @. canvas_color[channel, :, :] = canvas * color[channel, :, :]
        end

        canvas = canvas_color
    end

    return canvas
end

function get_blur_kernel(p::Float64)::Tuple{Vector{Float64},Int}
    i_max = Int(floor(p + 0.5))
    kernel = zeros(2 * i_max + 1)

    i = -i_max:i_max
    arg_right = @. (i[2:end-1] + 0.5) / p
    arg_left = @. (i[2:end-1] - 0.5) / p

    kernel[2:end-1] = @. 15 / 16 * (
        1 / 5 * (arg_right^5 - arg_left^5) +
        -2 / 3 * (arg_right^3 - arg_left^3) +
        arg_right - arg_left
    )

    x = (i_max - 0.5) / p

    c_outter = 0.5 - 15 / 16 * (x^5 / 5 - 2 * x^3 / 3 + x)

    kernel[1] = c_outter
    kernel[end] = c_outter

    return kernel, i_max
end

function add_depth_of_field(
    canvas::Array{Float64,3},
    t_int::Array{Float64,2},
    focus_curve,
)::Array{Float64,3}

    _, w, h = size(canvas)
    canvas_new = zero(canvas)

    for i = 1:w
        for j = 1:h
            if isinf(t_int[i, j])
                continue
            end
            focus = focus_curve(t_int[i, j])
            kernel, i_max = get_blur_kernel(focus)
            for i_ = -i_max:i_max
                for j_ = -i_max:i_max
                    i_abs = i + i_
                    j_abs = j + j_
                    if (i_abs < 1) || (i_abs > w) || (j_abs < 1) || (j_abs > h)
                        continue
                    end
                    @. canvas_new[:, i_abs, j_abs] +=
                        canvas[:, i, j] * kernel[i_max+i_+1] * kernel[i_max+j_+1]
                end
            end
        end
    end

    return canvas_new
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
    screen_size = [0.2, 0.1]
    screen_dist = [0.2]
    screen_res = [1800, 900]
    dropoff_curve(t) = clamp(1 - 0.3 * (t - 1), 0, 1)
    focus_curve(t) = 0.5 + 20 * abs(t - 3)
    function warp!(v::Vector{Float64})::Nothing
        v[3] = v[3] + 0.1 * sin(250 * v[2])
        v[1] = v[1] + 0.1 * sin(250 * v[2])
        return nothing
    end
    cam = Camera(
        [loc],
        [dir],
        [up],
        [screen_size],
        [screen_dist],
        [screen_res];
        warp = [warp!],
    )
    look_at!(cam, 1, [2.0, 2.0, 2.0], zeros(3))

    sphere = Sphere(zeros(3), 1.0)
    cube = Cube(zeros(3), 0.5)
    sponge = Menger_sponge(zeros(3), 0.5, 4)

    collect_metadata = Dict(:dim_int => Int)
    t_int, metadata = shape_view(cam, 1, sponge; collect_metadata)
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

    canvas = cam_is_source(t_int, dropoff_curve; color)
    canvas = add_depth_of_field(canvas, t_int, focus_curve)
    canvas = permutedims(canvas, [1, 3, 2])
    reverse!(canvas)

    colorview(RGB, canvas)
end

# Making animation
function warp_animation()
    julia_green = [0.22, 0.596, 0.149]
    julia_purple = [0.584, 0.345, 0.698]
    julia_red = [0.796, 0.235, 0.2]

    julia_colors = [julia_green, julia_purple, julia_red]

    loc = zeros(3)
    dir = zeros(3)
    up = zeros(3)
    screen_size = [0.1, 0.1]
    screen_dist = [0.2]
    screen_res = [1000, 1000]
    dropoff_curve(t) = clamp(1 - 0.3 * (t - 1), 0, 1)
    focus_curve(t) = 0.5 + 20 * abs(t - 3)
    function warp!(v::Vector{Float64})::Nothing
        v[3] = v[3] + 0.1 * sin(250 * v[2])
        v[1] = v[1] + 0.1 * sin(250 * v[2])
        return nothing
    end
    cam = Camera(
        [loc],
        [dir],
        [up],
        [screen_size],
        [screen_dist],
        [screen_res];
        warp = [warp!],
    )
    sponge = Menger_sponge(zeros(3), 0.5, 4)

    n_frames = 100
    framerate = 25
    ϕ = π / 3
    R = 4

    encoder_options = (; crf = 0)

    open_video_out(
        "C:/Users/bart1/Documents/Julia projects/Rays/Menger_sponge_spin.mp4",
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

            look_at!(cam, 1, from, to)

            collect_metadata = Dict(:dim_int => Int)
            t_int, metadata = shape_view(cam, 1, sponge; collect_metadata)
            color = zeros(3, screen_res...)
            for channel = 1:3
                color_channel = view(color, channel, :, :)
                @. color_channel[metadata[:dim_int]==channel] = 1.0
            end

            canvas = cam_is_source(t_int, dropoff_curve; color)
            canvas = permutedims(canvas, [1, 3, 2])
            reverse!(canvas)
            canvas = RGB{N0f8}.([canvas[channel, :, :] for channel = 1:3]...)
            write(writer, colorview(RGB, canvas))
        end
    end
end

simple_view()
