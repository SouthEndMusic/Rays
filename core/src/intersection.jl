abstract type Intersection end

struct Intersection_plain{F<:AbstractFloat} <: Intersection
    t::Vector{F}
end

struct Intersection_with_dim{F<:AbstractFloat,I<:Integer} <: Intersection
    t::Vector{F}
    dim::Vector{I}
end

struct Intersection_with_face{F<:AbstractFloat,I<:Integer} <: Intersection
    t::Vector{F}
    face::Vector{I}
end

get_default_intersection(::Sphere{F}) where {F} = Intersection_plain{F}([Inf])
get_default_intersection(::Cube{F}) where {F} = Intersection_with_dim{F,Int64}([Inf], [0])
get_default_intersection(shape::FractalShape) = get_default_intersection(first(shape.subshapes))
get_default_intersection(::TriangleShape{F,I}) where {F,I} =
    Intersection_with_face{F,I}([Inf], [0])

function set_default_intersection!(intersection::Intersection_plain)::Nothing
    intersection.t[1] = Inf
    return nothing
end

function set_default_intersection!(intersection::Intersection_with_dim)::Nothing
    intersection.t[1] = Inf
    intersection.dim[1] = 0
    return nothing
end

function set_default_intersection!(intersection::Intersection_with_face)::Nothing
    intersection.t[1] = Inf
    intersection.face[1] = 0
    return nothing
end

"""
Compute the intersection of a ray with a sphere as the smallest
real solution to a quadratic polynomial, if it exists.
"""
function intersect!(
    ray::Ray{F},
    sphere::Sphere;
    intersection::Union{Intersection_plain{F}, Nothing} = nothing,
)::Tuple{Bool,Intersection_plain{F}} where {F<:AbstractFloat}
    (; loc, dir) = ray
    (; center, Rsq) = sphere

    if isnothing(intersection)
        intersection = get_default_intersection(sphere)
    end

    diff = loc - center
    a = norm(dir)^2
    b = 2 * dot(dir, diff)
    c = norm(diff)^2 - Rsq
    discr = b^2 - 4 * a * c

    closer_intersection_found = false

    if discr >= 0
        t_int_candidate = (-b - sqrt(discr)) / (2 * a) 
        if t_int_candidate < intersection.t[1]
            closer_intersection_found = true
            intersection.t[1] = t_int_candidate
        end
    end
    
    return closer_intersection_found, intersection
end

"""
Compute the intersection of a ray with a cube by computing the intersections
with each of the 6 face planes and then checking whether the intersection is within the face.
"""
function intersect!(
    ray::Ray{F},
    cube::Cube{F};
    intersection::Union{Intersection_with_dim{F,<:Integer}, Nothing} = nothing,
)::Tuple{Bool,Intersection_with_dim{F,<:Integer}} where {F<:AbstractFloat}

    if isnothing(intersection)
        intersection = get_default_intersection(cube)
    end

    closer_intersection_found = false

    for dim = 1:3
        bound_small = cube.center[dim] - cube.R
        diff_bound_small = bound_small - ray.loc[dim]
        dir_dim_positive = (ray.dir[dim] > 0) # dir_dim = 0 not taken into account

        if diff_bound_small > 0.0
            if dir_dim_positive
                t_int_candidate = diff_bound_small / ray.dir[dim]
            else
                return intersection
            end
        else
            bound_big = cube.center[dim] + cube.R
            diff_bound_big = bound_big - ray.loc[dim]

            if diff_bound_big > 0.0
                if dir_dim_positive
                    t_int_candidate = diff_bound_big / ray.dir[dim]
                else
                    t_int_candidate = -diff_bound_small / ray.dir[dim]
                end
            else
                if dir_dim_positive
                    return intersection
                else
                    t_int_candidate = diff_bound_big / ray.dir[dim]
                end
            end
        end

        if t_int_candidate < intersection.t[1]
            other_dim = 0
            candidate = true

            while candidate && other_dim < 3
                other_dim += 1

                if other_dim !== dim
                    loc_int_other_dim_1 = ray.loc[other_dim] + t_int_candidate * ray.dir[other_dim]
                    if loc_int_other_dim_1 > cube.center[other_dim] + cube.R
                        candidate = false
                        continue
                    elseif loc_int_other_dim_1 < cube.center[other_dim] - cube.R
                        candidate = false
                        continue
                    end
                end
            end

            if candidate
                closer_intersection_found = true
                intersection.t[1] = t_int_candidate
                intersection.dim[1] = dim
            end
        end
    end
    return closer_intersection_found, intersection
end


"""
Compute the intersection of a ray with a fractal shape.
This is done recursively until the recursion depth of the fractal shape.
To compute the intersection of a ray with a subshape, the ray location is transformed.
"""
function intersect!(
    ray::Ray{F},
    fractal_shape::FractalShape{F};
    intersection::Union{Intersection, Nothing} = nothing,
    current_depth::Int = 0,
)::Tuple{Bool,Intersection} where {F <: AbstractFloat}
    (; subshapes, depth, shrink_factor) = fractal_shape

    if isnothing(intersection)
        intersection = get_default_intersection(first(fractal_shape.subshapes))
    end

    subshape_intersect = nothing

    for subshape in subshapes
        closer_intersection_found = intersect!(ray, subshape; intersection)[1]

        if closer_intersection_found
            subshape_intersect = subshape
        end
    end

    if !isnothing(subshape_intersect) && current_depth < depth
        set_default_intersection!(intersection)
        loc_transformed =
            shrink_factor * (ray.loc + fractal_shape.center - subshape_intersect.center)
        ray_transformed = Ray(loc_transformed, ray.dir)
        intersection.t[1] *= shrink_factor
        closer_intersection_found = intersect!(ray_transformed, fractal_shape; intersection, current_depth = current_depth + 1)[1]
        intersection.t[1] /= shrink_factor
        # intersection_candidate =
        #     intersect!(ray_transformed, fractal_shape; intersection, current_depth = current_depth + 1)
        # intersection_candidate.t[1] /= shrink_factor

        # if intersection_candidate.t < intersection.t
        #     intersection = intersection_candidate
        # end
    else
        closer_intersection_found = false
    end

    return closer_intersection_found, intersection
end

"""
Compute the intersection of a ray with a triangle given by
the triangle vertices.
"""
function intersect(
    ray::Ray{F},
    triangle_vertices::AbstractArray{F};
    normal::Union{AbstractVector{F}, Nothing} = nothing,
    t_int_prev::Union{F, Nothing} = nothing,
)::F where {F<:AbstractFloat}
    (; loc, dir) = ray

    if isnothing(normal)
        c = view(triangle_vertices, 3, :)
        u = view(triangle_vertices, 1, :) - c
        v = view(triangle_vertices, 2, :) - c
        normal = cross(u,v)
        normalize!(normal)
    end

    diff = view(triangle_vertices, 3, :) - loc
    t_int_candidate = dot(diff, normal)/dot(dir, normal)

    if t_int_candidate < 0.0
        return t_int_prev
    end

    if isnothing(t_int_prev)
        t_int_prev = convert(F, Inf)
    end

    if t_int_candidate ≥ t_int_prev
        return t_int_prev
    end

    if !isnothing(normal)
        c = view(triangle_vertices, 3, :)
        u = view(triangle_vertices, 1, :) - c
        v = view(triangle_vertices, 2, :) - c
    end

    u_normsq = dot(u, u)
    v_normsq = dot(v, v)
    inner_uv = dot(u, v)
    det = u_normsq * v_normsq - inner_uv^2
    diff2 = @. t_int_candidate * ray.dir - diff
    inner_u = dot(diff2, u)
    inner_v = dot(diff2, v)
    λ_1 = (inner_u * v_normsq - inner_v * inner_uv)/det

    if !(0 ≤ λ_1 ≤ 1)
        return t_int_prev
    end

    λ_2 = (inner_v * u_normsq - inner_u * inner_uv)/det

    if !(0 ≤ λ_2 ≤ 1)
        return t_int_prev
    end

    if λ_1 + λ_2 ≤ 1.0
        return t_int_candidate
    else
        return t_int_prev
    end
end

"""
Compute the intersection of a ray with a triangle shape
as the smallest intersection time over all triangles.
"""
function intersect!(
    ray::Ray{F},
    shape::TriangleShape{F};
    intersection::Union{Intersection_with_face{F,<:Integer},Nothing} = nothing,
)::Tuple{Bool,Intersection_with_face{F,<:Integer}} where {F<:AbstractFloat}
    (; vertices, faces, normals) = shape

    if isnothing(intersection)
        intersection = default_intersection(shape)
    end

    closer_intersection_found = false

    for i = 1:shape.n_faces
        triangle_vertices = view(vertices, view(faces, i, :), :)
        normal = view(normals, i, :)
        t_int_candidate = intersect(ray, triangle_vertices; normal, t_int_prev = intersection.t[1])
        if t_int_candidate < intersection.t[1]
            closer_intersection_found = true
            intersection.t[1] = t_int_candidate
            intersection.face[1] = i
        end
    end

    return closer_intersection_found, intersection
end
