"""
Create a Matrix{Float64} with intersection times of a shape per pixel.
If there is no intersection the intersection time is Inf.
Depending on the type of object, metadata of the intersections can be collected,
for instance the intersection dimension for cubes.
"""
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
                # TODO: Create a custom error for when the required metadata does not exist
                metadata[s][i, j] = getfield(int_metadata, s)
            end
        end
    end

    return t_int, metadata
end

"""
First create a grayscale canvas of the intersection times using the dropoff_curve;
this curve which maps positive numbers to [0,1] determines the brightness of a pixel
based on its distance to the camera. 
If a color Array{Float64,3} is provided, this is multiplied by the grayscale canvas to
produce a colored image with varying brightness.
"""
function cam_is_source(
    intersections::Matrix{Float64},
    dropoff_curve;
    color::Union{Array{Float64,3},Nothing} = nothing,
)::Array{Float64}

    where_intersect = .!isinf.(intersections)
    t_intersect = intersections[where_intersect]

    canvas = zeros(Float64, size(intersections)...)
    canvas[where_intersect] = @. dropoff_curve(t_intersect)

    # TODO: Move this to separate function for coloring
    if !isnothing(color)
        canvas_color = zeros(Float64, 3, size(intersections)...)

        for channel = 1:3
            @. canvas_color[channel, :, :] = canvas * color[channel, :, :]
        end

        canvas = canvas_color
    end

    return canvas
end

"""
Create a discrete distribution where the spread is based on p.
The distribution is based on integrating over pixels [i-1/2,i+1/2] ∩ [-p,p]
the continuous distribution 
f(x) = 15/16 (x^2-1)^2, x ∈ [-p,p].
"""
function get_blur_kernel(p::Float64)::Tuple{Vector{Float64},Int}
    Δi_max = Int(floor(p + 0.5))
    kernel = zeros(2 * Δi_max + 1)

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
            kernel, Δi_max = get_blur_kernel(focus)
            for Δi = -Δi_max:Δi_max
                for Δj = -Δi_max:Δi_max
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

    return canvas_new
end
