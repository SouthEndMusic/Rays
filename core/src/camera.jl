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
struct Camera{F<:AbstractFloat}
    loc::Vector{F}
    dir::Vector{F}
    up::Vector{F}
    right::Vector{F}
    screen_size::Vector{F} # height, width
    screen_dist::Vector{F}
    screen_res::Vector{Int} # height, width
    warp::Function
    function Camera(
        loc::Vector{F},
        dir::Vector{F},
        up::Vector{F},
        right::Vector{F},
        screen_size::Vector{F},
        screen_dist::Vector{F},
        screen_res,
        warp,
    ) where {F}
        return new{F}(loc, dir, up, right, screen_size, screen_dist, screen_res, warp)
    end
end

"""
Construct a camera object where the 'right' vector is computed using the cross product.
"""
function Camera(loc, dir, up, screen_size, screen_dist, screen_res; warp = nothing)::Camera
    right = cross(dir, up)

    warp = isnothing(warp) ? identity : warp

    return Camera(loc, dir, up, right, screen_size, screen_dist, screen_res, warp)
end

"""
Get a default camera instance
"""
function Camera(; float_type = Float32)::Camera{<:AbstractFloat}
    default_values_float =
        Vector{float_type}[zeros(3), [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.1, 0.1], [0.1]]
    screen_res = [100, 100]

    return Camera(default_values_float..., screen_res; warp = identity)
end


"""
Point the camera from 'from' to 'to'. 
'up' is defined to be a linear combination of unit length of
e_z = [0,0,1] and di,r and orthogonal to dir. This does not work if dir and e_z are proportional,
e.g. the camera points straight up or straight down.
"""
function look_at!(
    camera::Camera{F},
    from::AbstractVector{F},
    to::AbstractVector{F},
)::Nothing where {F}
    (; loc, dir, up, right) = camera

    loc .= from
    dir .= normalize(to - from)

    dir_z = dir[3]
    denom = sqrt(1 - dir_z^2)

    if denom ≈ 0.0
        error("In 'look_at!', the camera cannot point straight up or down.")
    end

    up .= -dir_z * dir
    up[3] += 1.0
    up ./= denom
    right .= cross(dir, up)

    return nothing
end

"""
Point the camera to 'to' from 'to' plus the spherical coordinates given by dist, θ and ϕ.
"""
function look_at!(
    camera::Camera{F},
    to::AbstractVector{F},
    dist::AbstractFloat,
    θ::AbstractFloat,
    ϕ::AbstractFloat,
)::Nothing where {F}
    from = to + dist * [
        cos(θ) * sin(ϕ)
        sin(θ) * sin(ϕ)
        cos(ϕ)
    ]
    from = convert(Vector{F}, from)
    look_at!(camera, from, to)
    return nothing
end

"""
Get an Array{Float64} with the resolution of the provided camera.
"""
function get_canvas(camera::Camera{F}; color::Bool = false)::Array{F} where {F}

    return if color
        zeros(F, 3, camera.screen_res...)
    else
        zeros(F, camera.screen_res...)
    end
end

"""
loc: The origin of the ray 
dir: The direction of the ray (unit vector)
"""
struct Ray{F<:AbstractFloat}
    loc::Vector{F}
    dir::Vector{F}
end

"""
Construct a ray.
"""
function Ray()::Ray
    return Ray(zeros(3), [1.0, 0.0, 0.0])
end

"""
Convert between rays with different type parameters.
"""
function Base.convert(::Type{Ray{F}}, ray::Ray) where {F}
    return Ray{F}(ray.loc, ray.dir)
end

"""
Compute the location and direction of a ray emited from the camera for 
the given pixel indices.
"""
function set_ray!(
    ray::Ray{F},
    camera::Camera{F},
    pixel_indices::Tuple{Int,Int},
)::Nothing where {F}
    (; screen_dist, screen_size, screen_res, dir, loc, up, right, warp) = camera

    screen_dist = screen_dist[1]

    if !(all(pixel_indices .> 0) && all(pixel_indices .<= screen_res))
        error("Pixel indices must fall inside screen res $screen_res, got $pixel_indices.")
    end

    i, j = pixel_indices
    s_h, s_w = screen_size
    res_h, res_w = screen_res

    @. ray.loc = loc + screen_dist * dir
    @. ray.loc -= ((i - 1) / (res_h - 1) - 0.5) * s_h * up
    @. ray.loc += ((j - 1) / (res_w - 1) - 0.5) * s_w * right

    @. ray.dir = ray.loc - loc
    normalize!(ray.dir)
    warp(ray.loc)

    return nothing
end
