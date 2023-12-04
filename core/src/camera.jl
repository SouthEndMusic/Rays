"""
name: the name of the camera
loc: The location of the camera in space
dir: The direction the camera is pointing in (vector of unit length)
up: The upwards direction 
right: The right direction
screen_size: The size of the screen in world units
screen_dist: The distance between loc and the image plane
screen_res: The resolution of the resulting render
t_intersect: The intersection time per pixel
canvas: The array to which the render is written
warp!: A function which changes the origin of rays as they leave the camera
dropoff_curve: Function [0,∞) → [0,1] for the brightness dropoff with the distance to the camera
	when using cam_is_source!
focus_curve: Function [0,∞) → R for the blurriness as function of the intersection t

dir, up and right are orthormal and completely fix the orientation of the camera.
"""
struct Camera{F<:AbstractFloat,FC<:Union{ScalarFunc{F},Nothing},MF<:AbstractMatrix{F}}
    name::Symbol
    loc::Vector{F}
    dir::Vector{F}
    up::Vector{F}
    right::Vector{F}
    screen_size::Vector{F} # height, width
    screen_dist::Vector{F}
    t_intersect::Matrix{F}
    canvas::Array{F,3}
    warp!::Transform{F,MF}
    dropoff_curve::ScalarFunc{F}
    focus_curve::FC
end

function Base.show(io::IO, camera::Camera)::Nothing
    (; name) = camera
    print(io, "<Camera '$name'>")
    return nothing
end

function get_screen_res(camera::Camera)::Tuple{Int,Int}
    return size(camera.t_intersect)
end

"""
Set the focus curve of a camera.
Note: this is not a mutating function, it creates a new Camera instance.
"""
function set_focus_curve(
    camera::Camera{F},
    focus_curve::Union{Function,Nothing},
)::Camera{F} where {F}
    if !isnothing(focus_curve)
        focus_curve = ScalarFunc{F}(focus_curve)
    end
    return @set camera.focus_curve = focus_curve
end

"""
Set the warp! of a camera.
Note: this is not a mutating function, it creates a new Camera instance.
"""
function set_warp(camera::Camera{F,FC,MF}, warp!::Function)::Camera{F} where {F,FC,MF}
    warp! = Transform{F,MF}(warp!)
    return @set camera.warp! = warp!
end

"""
Construct a camera object where the 'right' vector is computed using the cross product.
"""
function Camera(
    loc::VF,
    dir::VF,
    up::VF,
    screen_size::VF,
    screen_dist::VF,
    screen_res::Tuple{Int,Int};
    warp!::Union{Function,Nothing} = nothing,
    name::Union{Symbol,Nothing} = nothing,
    dropoff_curve::Union{Function,Nothing} = nothing,
    focus_curve::Union{Function,Nothing} = nothing,
)::Camera where {VF<:AbstractVector{F} where {F<:AbstractFloat}}
    F = eltype(VF)
    right = cross(dir, up)

    warp! = isnothing(warp!) ? identity : warp!

    canvas = similar(loc, (3, screen_res...))
    t_intersect = similar(loc, screen_res)

    # Generate default name if not given
    if isnothing(name)
        name = snake_case_name(Camera)
    end

    if isnothing(dropoff_curve)
        dropoff_curve = ScalarFunc{F}(t -> convert(F, 1.0))
    end

    if !isnothing(focus_curve)
        focus_curve = ScalarFunc{F}(focus_curve)
    end

    MF = typeof(similar(loc, (3, 3)))

    return Camera(
        name,
        loc,
        dir,
        up,
        right,
        screen_size,
        screen_dist,
        t_intersect,
        canvas,
        Transform{F,MF}(warp!),
        ScalarFunc{F}(dropoff_curve),
        focus_curve,
    )
end

"""
Get a default camera instance
"""
function Camera(;
    screen_res::Tuple{Int,Int} = (100, 100),
    screen_size::Vector{<:AbstractFloat} = [0.1, 0.1],
    vector_prototype::AbstractVector{F} where {F<:AbstractFloat} = zeros(Float32, 3),
    name::Union{Symbol,Nothing} = nothing,
)::Camera
    default_values_float = typeof(vector_prototype)[
        zeros(3),
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        screen_size,
        [0.1],
    ]

    return Camera(default_values_float..., screen_res; name)
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

    copyto!(loc, from)
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
    dist::F,
    θ::F,
    ϕ::F,
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
loc: The origin of the ray 
dir: The direction of the ray (unit vector)
"""
struct Ray{MF<:AbstractMatrix{F} where {F<:AbstractFloat}}
    loc::MF
    dir::MF
end

"""
Construct a ray.
"""
function Ray(
    n::Int;
    matrix_prototype::AbstractMatrix{F} where {F<:AbstractFloat} = zeros(Float32, 3, 3),
)::Ray
    return Ray(similar(matrix_prototype, (n, 3)), similar(matrix_prototype, (n, 3)))
end

"""
Compute the location and direction of a ray emited from the camera for 
the given pixel indices.
"""
function set_ray!(
    ray_loc::AbstractVector{F},
    ray_dir::AbstractVector{F},
    camera::Camera{F},
    pixel_indices::Tuple{Int,Int},
)::Nothing where {F}
    (; screen_dist, screen_size, dir, loc, up, right, warp!) = camera
    screen_res = get_screen_res(camera)

    screen_dist = screen_dist[1]

    i, j = pixel_indices
    s_h, s_w = screen_size
    res_h, res_w = screen_res

    @. ray_loc = loc + screen_dist * dir
    @. ray_loc -= ((i - 1) / (res_h - 1) - 0.5) * s_h * up
    @. ray_loc += ((j - 1) / (res_w - 1) - 0.5) * s_w * right

    @. ray_dir = ray_loc - loc
    normalize!(ray_dir)
    warp!(ray_loc)


    return nothing
end
