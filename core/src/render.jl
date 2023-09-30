"""
Create a Matrix{Float64} with intersection times of a shape per pixel.
If there is no intersection the intersection time is Inf.
Depending on the type of object, metadata of the intersections can be collected,
for instance the intersection dimension for cubes.
"""
function shape_view(
    camera::Camera{F},
    shape::Shape{F};
    data_variables::Vector{Symbol} = Symbol[],
)::NamedTuple where {F}

    (; screen_res) = camera
    if !(:t in data_variables)
        push!(data_variables, :t)
    end

    intersection_example::Intersection{F} = Intersection()
    data_types = DataType[]
    for data_var in data_variables
        if hasfield(Intersection, data_var)
            push!(data_types, eltype(getfield(intersection_example, data_var)))
        else
            error("Intersections cannot have metadata $data_var.")
        end
    end

    data_matrices = [zeros(T, screen_res...) for T in data_types]
    data = NamedTuple{Tuple(data_variables)}(Tuple(data_matrices))

    intersections = Intersection{F}[Intersection() for i = 1:nthreads()]

    @threads for I in CartesianIndices(data.t)
        intersection = intersections[threadid()]
        set_ray!(intersection.ray, camera, Tuple(I))
        intersect!(intersection, shape)

        for data_var in data_variables
            getfield(data, data_var)[I] = getfield(intersection, data_var)[1]
        end
        reset_intersection!(intersection)
    end

    return data
end

"""
First create a grayscale canvas of the intersection times using the dropoff_curve;
this curve which maps positive numbers to [0,1] determines the brightness of a pixel
based on its distance to the camera. 
If a color Array{Float64,3} is provided, this is multiplied by the grayscale canvas to
produce a colored image with varying brightness.
"""
function cam_is_source(
    t_int::AbstractMatrix{F};
    dropoff_curve::Union{Function,Nothing} = nothing,
)::AbstractMatrix{F} where {F}

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

    canvas_grayscale = zeros(F, size(t_int)...)
    canvas_grayscale[where_intersect] = @. curve(t_int_noninf)

    return canvas_grayscale
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
function add_depth_of_field(
    canvas::AbstractArray{F,3},
    t_int::AbstractMatrix{F},
    focus_curve,
)::AbstractArray{F,3} where {F}

    _, h, w = size(canvas)
    canvas_new = zeros(F, size(canvas)...)

    @threads for I in CartesianIndices(t_int)
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
                    canvas[:, i, j] * kernel[Δi_max+Δi+1] * kernel[Δi_max+Δj+1]
            end
        end
    end

    return canvas_new
end

"""
Apply a colormapping metadata -> colormapping to the color array.
color_palette: (3, n_colors) array
metadata: a metadata value of n yields the nth color in color color_palette
a metadata value of 0 yields no change in color
"""
function add_color!(
    color::AbstractArray{F,3},
    color_palette::AbstractMatrix{F},
    metadata::AbstractMatrix{Int},
)::Nothing where {F}
    @threads for I in CartesianIndices(metadata)
        metadata_value = metadata[I]
        if !iszero(metadata_value)
            color[:, I] = view(color_palette, :, metadata_value)
        end
    end
    return nothing
end

"""
Multiply a grayscale canvas by a color array to get a color canvas.
"""
function apply_color(
    canvas_grayscale::AbstractMatrix{F},
    color::AbstractArray{F},
)::AbstractArray{F,3} where {F}
    res_canvas_grayscale = size(canvas_grayscale)
    res_color = size(color)[2:3]
    @assert res_canvas_grayscale == res_color "canvas_grayscale and color must have the same resolution, got $res_canvas_grayscale and $res_color respectively."

    canvas_color = zeros(F, size(color)...)

    @threads for I in CartesianIndices(canvas_grayscale)
        @. canvas_color[:, I] = canvas_grayscale[I] * color[:, I]
    end

    return canvas_color
end
