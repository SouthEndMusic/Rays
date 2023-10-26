"""
Create a Matrix{Float64} with intersection times of a shape per pixel.
If there is no intersection the intersection time is Inf.
Depending on the type of object, metadata of the intersections can be collected,
for instance the intersection dimension for cubes.
"""
function shape_view!(
    scene::Scene{F};
    camera_index::Int = 1,
    shape_index::Int = 1,
)::Nothing where {F}

    camera = scene.cameras[camera_index]
    (; intersection_data_float, intersection_data_int) = camera
    data_variables_float = keys(intersection_data_float)
    data_variables_int = keys(intersection_data_int)

    shape = scene.shapes[shape_index]
    intersections = Intersection{F}[Intersection() for i = 1:nthreads()]

    n_chunks = 10 * nthreads()
    CI = CartesianIndices(Tuple(camera.screen_res))
    numel = prod(camera.screen_res)

    @threads for c in 1:n_chunks
        for I_flat in c:n_chunks:numel
            I = CI[I_flat]
            intersection = intersections[threadid()]
            set_ray!(intersection.ray, camera, Tuple(I))
            intersect!(intersection, shape)

            # These are written out explicitly per variable
            # because looping leads to runtime-dispatch
            if :t in data_variables_float
                intersection_data_float[:t][I] = intersection.t[1]
            end
            if :dim in data_variables_int
                intersection_data_int[:dim][I] = intersection.dim[1]
            end
            if :face in data_variables_int
                intersection_data_int[:face][I] = intersection.face[1]
            end
            reset_intersection!(intersection)
        end
    end

    return nothing
end

"""
First create a grayscale canvas of the intersection times using the dropoff_curve;
this curve which maps positive numbers to [0,1] determines the brightness of a pixel
based on its distance to the camera. 
If a color Array{Float64,3} is provided, this is multiplied by the grayscale canvas to
produce a colored image with varying brightness.
"""
function cam_is_source!(
    cam::Camera{F};
    dropoff_curve::Union{Function,Nothing} = nothing,
)::Nothing where {F}
    (; canvas_grayscale, intersection_data_float) = cam
    t_int = intersection_data_float[:t]

    # Create a dropoff curve from the closest intersection to the
    # furthest intersection fi no dropoff_curve is given
    if isnothing(dropoff_curve)
        Max = maximum(x -> isinf(x) ? -Inf : x, t_int)
        Min = minimum(t_int)
        Diff = Max - Min
        curve(x) = 1 - (x - Min) / Diff
    else
        curve = dropoff_curve
    end

    where_intersect = .!isinf.(t_int)
    t_int_noninf = t_int[where_intersect]

    canvas_grayscale .= 0.0
    canvas_grayscale[where_intersect] = @. curve(t_int_noninf)

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
function add_depth_of_field!(camera::Camera{F}, focus_curve::Function)::Nothing where {F}

    (; canvas_color, intersection_data_float) = camera

    _, h, w = size(canvas_color)
    t_int = intersection_data_float[:t]
    canvas_new = zeros(F, size(canvas_color)...)

    n_chunks = 10 * nthreads()
    CI = CartesianIndices(Tuple(camera.screen_res))
    numel = prod(camera.screen_res)

    @threads for c in 1:n_chunks
        for I_flat in c:n_chunks:numel
            I = CI[I_flat]
            if isinf(t_int[I])
                continue
            end
            focus = focus_curve(t_int[I])
            kernel, Δi_max = get_blur_kernel(focus)
            i, j = Tuple(I)
            for Δi = -Δi_max:Δi_max
                for Δj = -Δi_max:Δi_max
                    i_abs = i + Δi
                    j_abs = j + Δj
                    if (i_abs < 1) || (i_abs > w) || (j_abs < 1) || (j_abs > h)
                        continue
                    end
                    @. canvas_new[:, i_abs, j_abs] +=
                        canvas_color[:, i, j] * kernel[Δi_max+Δi+1] * kernel[Δi_max+Δj+1]
                end
            end
        end
    end

    canvas_color .= canvas_new
    return nothing
end

"""
Apply a colormapping metadata -> colormapping to the color array.
color_palette: (3, n_colors) array
metadata: a metadata value of n yields the nth color in color color_palette
a metadata value of 0 yields no change in color
"""
function get_color!(
    camera::Camera{F},
    variable::Symbol,
    color_palette::AbstractMatrix{F},
)::Nothing where {F}
    (; color, intersection_data_int) = camera
    if variable ∉ keys(intersection_data_int)
        error("Integer intersection data of variable '$variable' has not been collected.")
    end
    data = intersection_data_int[variable]
    n_chunks = 10 * nthreads()
    CI = CartesianIndices(Tuple(camera.screen_res))
    numel = prod(camera.screen_res)
    @threads for c in 1:n_chunks
        for I_flat in c:n_chunks:numel
            I = CI[I_flat]
            data_value = data[I]
            if !iszero(data_value)
                color[:, I] = view(color_palette, :, data_value)
            end
        end
    end
    return nothing
end

"""
Multiply a grayscale canvas by a color array to get a color canvas.
"""
function apply_color!(camera::Camera{F})::Nothing where {F}
    (; canvas_grayscale, color, canvas_color) = camera
    n_chunks = 10 * nthreads()
    CI = CartesianIndices(Tuple(camera.screen_res))
    numel = prod(camera.screen_res)
    @threads for c in 1:n_chunks
        for I_flat in c:n_chunks:numel
            I = CI[I_flat]
            canvas_color[:, I] .= canvas_grayscale[I]
            @views(canvas_color[:, I] .*= color[:, I])
    
        end
    end
    return nothing
end
