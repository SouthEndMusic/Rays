function set_dropoff_curve_default!(scene::Scene{F}, camera::Camera{F})::Camera{F} where {F}
    if length(scene.shapes) == 0
        error("Cannot determine default dropoff curve without shapes in the scene.")
    end
    dist_max = zero(F)
    for shape in values(scene.shapes)
        dist_max = max(dist_max, norm(camera.loc - shape.center))
    end
    dropoff_curve =
        ScalarFunc{F}(t -> max(zero(F), one(F) - t / (convert(F, 1.5) * dist_max)))

    camera_new = @set camera.dropoff_curve = dropoff_curve
    scene.cameras[camera.name] = camera_new
    return camera_new
end

function cam_is_source!(
    camera::Camera{F},
    t_intersect::F,
    name_intersected::Symbol,
    pixel_indices::Tuple{Int,Int},
)::Nothing where {F}
    (; canvas, dropoff_curve) = camera

    if name_intersected == :none
        shade = zero(F)
    else
        shade = dropoff_curve(t_intersect)
    end

    canvas[:, pixel_indices...] .= shade
    return nothing
end

"""
Create a discrete distribution where the spread is based on p.
The distribution is based on integrating over pixels [i-1/2,i+1/2] ∩ [-p,p]
the continuous distribution 
f(x) = 15/16 (x^2-1)^2, x ∈ [-p,p].
"""
function get_blur_kernel(p::F)::Tuple{AbstractVector{F},Int} where {F}
    Δi_max = Int(floor(p + 0.5))
    kernel = zeros(F, 2 * Δi_max + 1)

    Δi = -Δi_max:Δi_max
    arg_right = @. (Δi[2:end-1] + 0.5) / p
    arg_left = @. (Δi[2:end-1] - 0.5) / p

    kernel[2:end-1] = @. 15 / 16 * (
        1 / 5 * (arg_right^5 - arg_left^5) +
        -2 / 3 * (arg_right^3 - arg_left^3) +
        arg_right - arg_left
    )

    x = (Δi_max - 0.5) / p
    c_outter = 0.5 - 15 / 16 * (x^5 / 5 - 2 * x^3 / 3 + x)

    kernel[1] = c_outter
    kernel[end] = c_outter

    return kernel, Δi_max
end

"""
The depth of field effect is created by applying a kernel to pixel to 
'smear out' its value over the neighbouring pixels, where the spread
is given by the focus curve at the intersection time at that pixel.
"""
function add_depth_of_field!(camera::Camera{F})::Nothing where {F}

    (; canvas, t_intersect, focus_curve) = camera

    _, h, w = size(canvas)
    canvas_new = zeros(F, size(canvas)...)

    n_chunks = 10 * nthreads()
    CI = CartesianIndices(Tuple(camera.screen_res))
    numel = prod(camera.screen_res)

    @threads for c ∈ 1:n_chunks
        for I_flat ∈ c:n_chunks:numel
            I = CI[I_flat]
            if isinf(t_intersect[I])
                continue
            end
            focus = focus_curve(t_intersect[I])
            kernel, Δi_max = get_blur_kernel(focus)
            i, j = Tuple(I)
            for Δi ∈ -Δi_max:Δi_max
                for Δj ∈ -Δi_max:Δi_max
                    i_abs = i + Δi
                    j_abs = j + Δj
                    if (i_abs < 1) || (i_abs > w) || (j_abs < 1) || (j_abs > h)
                        continue
                    end
                    @. canvas_new[:, i_abs, j_abs] +=
                        canvas[:, i, j] * kernel[Δi_max+Δi+1] * kernel[Δi_max+Δj+1]
                end
            end
        end
    end

    canvas .= canvas_new
    return nothing
end

"""
Apply a colormapping metadata -> colormapping to the color array.
color_palette: (3, n_colors) array
metadata: a metadata value of n yields the nth color in color color_palette
a metadata value of 0 yields no change in color
"""
function set_color!(
    camera::Camera{F},
    variable::Symbol,
    intersection::Intersection{F},
    color_palette::AbstractMatrix{F},
    pixel_indices::Tuple{Int,Int},
)::Nothing where {F}
    (; canvas) = camera

    data_value = getfield(intersection, variable)[1]
    name_intersected = intersection.name_intersected[1]

    if name_intersected !== :none
        @views(canvas[:, pixel_indices...] .*= color_palette[:, data_value])
    end
    return nothing
end

function render!(
    scene::Scene{F};
    name_camera::Union{Symbol,Nothing} = nothing,
    cam_is_source::Bool = true,
    color_palette::Union{AbstractMatrix{F},Nothing} = nothing,
    variable::Union{Symbol,Nothing} = nothing,
)::Nothing where {F}

    # If no camera is specified, take the first one
    if isnothing(name_camera)
        camera = first(values(scene.cameras))
    else
        camera = scene.cameras[name_camera]
    end

    # Number of tasks
    n_tasks = 10 * nthreads()

    # All indices of the render
    CI = CartesianIndices(Tuple(camera.screen_res))

    # Number of pixels
    n_pixels = prod(camera.screen_res)

    # Intersection objects
    intersections = Intersection{F}[Intersection() for i ∈ 1:nthreads()]

    # Reset the canvas
    camera.canvas .= one(F)

    @threads for task ∈ 1:n_tasks
        intersection = intersections[threadid()]
        for I_flat ∈ task:n_tasks:n_pixels
            indices = CI[I_flat]
            set_ray!(intersection.ray, camera, Tuple(indices))
            for shape in values(scene.shapes)
                intersect!(intersection, shape)
            end
            camera.t_intersect[indices] = only(intersection.t)
            if cam_is_source
                cam_is_source!(
                    camera,
                    only(intersection.t),
                    only(intersection.name_intersected),
                    Tuple(indices),
                )
            end
            if !isnothing(color_palette) && !isnothing(variable)
                set_color!(camera, variable, intersection, color_palette, indices.I)
            end
            reset_intersection!(intersection)
        end
    end

    if !isnothing(camera.focus_curve)
        add_depth_of_field!(camera)
    end

    return nothing
end
