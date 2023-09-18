"""
Struct for possibly multiple cameras. Per camera:
loc: The location of the camera in space
dir: The direction the camera is pointing in (vector of unit length)
up: The upwards direction 
right: The right direction
screen_size: The size of the screen in world units
screen_dist: The distance between loc and the image plane
screen_res: The resolution of the resulting render
warp: A function which changes the origin of rays as they leave the camera

dir, up and right are orthormal and completely fix the orientation of the camera.
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
Construct a camera object where the 'right' vectors are computed using the cross product
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
Get a trivial camera instance
"""
function Camera(len::Int)::Camera
    default_values =
        [zeros(3), [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.1, 0.1], [0.1], [100, 100]]
    return Camera([fill(default_value, len) for default_value in default_values]...)
end


"""
Point the camera with the given cam_index from 'from' to 'to'. 
'up' is defined to be a linear combination of unit length of
e_z = [0,0,1] and di,r and orthogonal to dir. This does not work if dir and e_z are proportional,
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

    if !(all(pixel_indices .> 0) && all(pixel_indices .<= screen_res))
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
