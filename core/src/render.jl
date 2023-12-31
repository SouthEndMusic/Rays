"""
Set the dropoff curve as:

f(t) = max(0, 1 - t/(1.5*dist_max)) 

where dist_max is the maximum distance between the camera location and
the center of a shape.
"""
function set_dropoff_curve_default!(scene::Scene{F}, camera::Camera{F})::Camera{F} where {F}
    if length(scene.shapes) == 0
        error("Cannot determine default dropoff curve without shapes in the scene.")
    end
    dist_max = zero(F)

    for transform ∈ values(scene.transforms)
        center = ismissing(transform.translation) ? zeros(F, 3) : transform.translation
        dist_max = max(dist_max, norm(camera.loc - center))
    end
    dropoff_curve =
        ScalarFunc{F}(t -> max(zero(F), one(F) - t / (convert(F, 1.5) * dist_max)))

    camera_new = @set camera.dropoff_curve = dropoff_curve
    scene.cameras[camera.name] = camera_new
    return camera_new
end

"""
Apply the dropoff curve.
"""
function cam_is_source!(
    camera::Camera{F},
    t::AbstractVector{F},
    pixel_indices::Tuple{Int,Int},
)::Nothing where {F}
    (; canvas, dropoff_curve) = camera
    t = t[1]

    if isinf(t)
        shade = zero(F)
    else
        shade = dropoff_curve(t)
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
    CI = CartesianIndices(get_screen_res(camera))
    numel = prod(get_screen_res(camera))

    @batch for c ∈ 1:n_chunks
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
Multiply the shade of the canvas by the color given by
the texture of the intersected shape.
If no intersection was found, set the color to black.
"""
function set_color!(
    color::AbstractVector{F},
    camera::Camera{F,FC,MF},
    ray_loc::AbstractVector{F},
    ray_dir::AbstractVector{F},
    t::AbstractVector{F},
    texturers::Dict,
    name_intersected::Symbol,
    pixel_indices::Tuple{Int,Int},
    cache_int::VI,
    cache_float::VF,
)::Nothing where {F,VF<:AbstractVector{F},VI<:AbstractVector{Int},FC,MF}
    (; canvas) = camera

    # Get texture color
    if name_intersected == :none
        color .= zero(F)
    else
        texturer! = texturers[name_intersected]
        texturer!(color, cache_int, cache_float, ray_loc, ray_dir, t)
    end

    # Set texture color
    @views(canvas[:, pixel_indices...] .*= color)
    return nothing
end

"""
Combine all intersectors into one function.
"""
function combine_intersectors(intersectors)::Function
    return intersectors!(
        t,
        ray_loc,
        ray_dir,
        ray_camera_loc,
        ray_camera_dir,
        cache_int,
        cache_float,
        name_intersected,
    ) = begin
        for intersector! ∈ values(intersectors)
            intersector!(
                t,
                ray_loc,
                ray_dir,
                ray_camera_loc,
                ray_camera_dir,
                cache_int,
                cache_float,
                name_intersected,
            )
        end
    end
end

"""
Intersect shapes efficiently based on the scene partition.
"""
function partition_intersector(scene::Scene)::Function
    (; partition, intersectors) = scene
    (; partition_nodes) = partition
    return intersector!(
        t,
        ray_loc,
        ray_dir,
        ray_camera_loc,
        ray_camera_dir,
        cache_int,
        cache_float,
        name_intersected,
    ) = begin
        nodes_to_process = [1]
        while !isempty(nodes_to_process)
            node_index = pop!(nodes_to_process)
            node = partition_nodes[node_index]
            if bounding_box_intersect(ray_camera_loc, ray_camera_dir, node.bounding_box)
                if isempty(node.child_indices)
                    # Doing this instead of looping over node.identifiers
                    # prevents runtime dispatch
                    for (name, intersector!) ∈ intersectors
                        if name ∉ node.identifiers
                            continue
                        end
                        intersector!(
                            t,
                            ray_loc,
                            ray_dir,
                            ray_camera_loc,
                            ray_camera_dir,
                            cache_int,
                            cache_float,
                            name_intersected,
                        )
                    end
                else
                    append!(nodes_to_process, node.child_indices)
                end
            end
        end
    end
end

"""
The main loop over pixels voor rendering.
"""
function render!(
    scene::Scene{F};
    name_camera::Union{Symbol,Nothing} = nothing,
    cam_is_source::Bool = true,
)::Nothing where {F}

    # If no camera is specified, take the first one
    camera::Camera{F} = if isnothing(name_camera)
        first(values(scene.cameras))
    else
        scene.cameras[name_camera]
    end

    (; intersectors, texturers, intersections) = scene
    (; canvas, t_intersect, focus_curve) = camera

    # Number of tasks
    n_tasks = 10 * nthreads()

    # All indices of the render
    CI = CartesianIndices(get_screen_res(camera))

    # Number of pixels
    n_pixels = prod(get_screen_res(camera))

    # Reset the canvas
    canvas .= one(F)

    # Combine all intersectors in one function, 
    # this prevents runtime dispatch
    if isempty(scene.partition.partition_nodes)
        intersectors! = combine_intersectors(intersectors)
    else
        intersectors! = partition_intersector(scene)
    end

    @batch for task ∈ 1:n_tasks
        ray_loc,
        ray_dir,
        ray_camera_loc,
        ray_camera_dir,
        t,
        cache_int,
        cache_float,
        color,
        name_intersected = get_caches(intersections, threadid())

        for I_flat ∈ task:n_tasks:n_pixels
            reset_intersection!(t, name_intersected)
            indices = CI[I_flat]

            # Set the ray based on the current pixel of the camera
            set_ray!(ray_camera_loc, ray_camera_dir, camera, Tuple(indices))

            # Calculate shape intersections
            intersectors!(
                t,
                ray_loc,
                ray_dir,
                ray_camera_loc,
                ray_camera_dir,
                cache_int,
                cache_float,
                name_intersected,
            )

            # Apply shading
            t_intersect[indices] = t[1]
            if cam_is_source
                cam_is_source!(camera, t, Tuple(indices))
            end

            set_color!(
                color,
                camera,
                ray_camera_loc,
                ray_camera_dir,
                t,
                texturers,
                name_intersected[1],
                Tuple(indices),
                cache_int,
                cache_float,
            )
        end
    end

    if !isnothing(focus_curve)
        add_depth_of_field!(camera)
    end

    return nothing
end
