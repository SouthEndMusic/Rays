struct Intersection{F<:AbstractFloat}
    ray::Ray{F}
    rays_transformed::Vector{Ray{F}}
    t::Vector{F}
    # For cube intersections
    dim::Vector{Int}
    # For triangle intersections
    u::Vector{F}
    v::Vector{F}
    diff::Vector{F}
    diff2::Vector{F}
    face::Vector{Int}
end

function Base.convert(::Type{Intersection{F}}, intersection::Intersection) where {F}
    return Intersection{F}(
        [getfield(intersection, fieldname) for fieldname in fieldnames(Intersection)]...,
    )
end

function Intersection()::Intersection
    return Intersection(
        Ray(), # ray
        Ray{Float64}[],
        [Inf], # t 
        [0], # dim 
        zeros(3), # u
        zeros(3), # v
        zeros(3), # diff
        zeros(3), # diff2
        [0], # face
    )
end

function reset_intersection!(intersection::Intersection)::Nothing
    intersection.t[1] = Inf
    return nothing
end

"""
Compute the intersection of a ray with a sphere as the smallest
real solution to a quadratic polynomial, if it exists.
"""
function intersect!(intersection::Intersection{F}, sphere::Sphere)::Bool where {F}
    (; loc, dir) = intersection.ray
    (; center, Rsq) = sphere

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

    return closer_intersection_found
end

"""
Compute the intersection of a ray with a cube by computing the intersections
with each of the 6 face planes and then checking whether the intersection is within the face.
"""
function intersect!(intersection::Intersection{F}, cube::Cube{F})::Bool where {F}
    (; ray) = intersection

    closer_intersection_found = false

    for dim ∈ 1:3
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
                    loc_int_other_dim_1 =
                        ray.loc[other_dim] + t_int_candidate * ray.dir[other_dim]
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
    return closer_intersection_found
end


"""
Compute the intersection of a ray with a fractal shape.
This is done recursively until the recursion depth of the fractal shape.
To compute the intersection of a ray with a subshape, the ray location is transformed.
"""
function intersect!(
    intersection::Intersection{F},
    fractal_shape::FractalShape{F};
    current_depth::Int = 1,
)::Bool where {F}
    (; ray) = intersection
    (; subshapes, depth, shrink_factor) = fractal_shape

    if current_depth == 1
        for i = 1:depth
            if length(intersection.rays_transformed) < i
                push!(intersection.rays_transformed, ray)
            end
        end
    end

    ray.loc .= intersection.rays_transformed[depth].loc
    subshape_intersect = nothing

    for subshape in subshapes
        closer_intersection_found = intersect!(intersection, subshape)

        if closer_intersection_found
            subshape_intersect = subshape
        end
    end

    if !isnothing(subshape_intersect) && current_depth < depth
        reset_intersection!(intersection)
        ray_transformed = intersection.rays_transformed[current_depth+1]
        ray_transformed.loc .= ray.loc
        ray_transformed.loc .+= fractal_shape.center
        ray_transformed.loc .-= subshape_intersect.center
        ray_transformed.loc .*= shrink_factor

        intersection.t[1] *= shrink_factor
        closer_intersection_found =
            intersect!(intersection, fractal_shape, current_depth = current_depth + 1)
        intersection.t[1] /= shrink_factor
    else
        closer_intersection_found = false
    end

    return closer_intersection_found
end

"""
Compute the intersection of a ray with a triangle given by
the triangle vertices.
"""
function intersect!(
    intersection::Intersection{F},
    triangle_vertices::AbstractArray{F};
    normal::Union{AbstractVector{F},Nothing} = nothing,
    t_int_prev::Union{F,Nothing} = nothing,
)::F where {F}
    (; ray, diff, diff2, u, v) = intersection
    (; loc, dir) = intersection.ray

    if isnothing(normal)
        c = view(triangle_vertices, 3, :)
        u .= view(triangle_vertices, 1, :) - c
        v .= view(triangle_vertices, 2, :) - c
        normal = cross(u, v)
        normalize!(normal)
    end

    diff .= 0.0
    diff .+= view(triangle_vertices, 3, :)
    diff .-= loc
    t_int_candidate = dot(diff, normal) / dot(dir, normal)

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
        u .= 0.0
        v .= 0.0
        u .+= view(triangle_vertices, 1, :)
        v .+= view(triangle_vertices, 2, :)
        u .-= c
        v .-= c
    end

    u_normsq = dot(u, u)
    v_normsq = dot(v, v)
    inner_uv = dot(u, v)
    det = u_normsq * v_normsq - inner_uv^2
    diff2 .= @. t_int_candidate * ray.dir - diff
    inner_u = dot(diff2, u)
    inner_v = dot(diff2, v)
    λ_1 = (inner_u * v_normsq - inner_v * inner_uv) / det

    if !(0 ≤ λ_1 ≤ 1)
        return t_int_prev
    end

    λ_2 = (inner_v * u_normsq - inner_u * inner_uv) / det

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
function intersect!(intersection::Intersection{F}, shape::TriangleShape{F};)::Bool where {F}
    (; vertices, faces, normals) = shape

    closer_intersection_found = false

    for i ∈ 1:shape.n_faces
        triangle_vertices = view(vertices, view(faces, i, :), :)
        normal = view(normals, i, :)
        t_int_candidate = intersect!(
            intersection,
            triangle_vertices;
            normal,
            t_int_prev = intersection.t[1],
        )
        if t_int_candidate < intersection.t[1]
            closer_intersection_found = true
            intersection.t[1] = t_int_candidate
            intersection.face[1] = i
        end
    end

    return closer_intersection_found
end
