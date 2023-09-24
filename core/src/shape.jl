abstract type Shape end
abstract type Intersection end

# Intersection types only contain one value per field, but they are vectors to make them mutable.

"""
Get default values for the metadata of the intersection of
a certain shape for when there is no intersection.
"""
struct Sphere{F<:AbstractFloat} <: Shape
    center::Vector{F}
    R::F
    Rsq::F
end

struct Intersection_plain{F<:AbstractFloat} <: Intersection
    t::Vector{F}
end

default_intersection(::Sphere{F}) where {F} = Intersection_plain{F}([Inf])
Sphere(center::Vector{F}, R::F) where {F<:AbstractFloat} = Sphere(center, R, R^2)

struct Cube{F<:AbstractFloat} <: Shape
    center::Vector{F}
    R::F
end

struct Intersection_with_dim{F<:AbstractFloat,I<:Integer} <: Intersection
    t::Vector{F}
    dim::Vector{I}
end

default_intersection(::Cube{F}) where {F} = Intersection_with_dim{F,Int64}([Inf], [0])

"""
A shape where each subshape is substituted by a shrinked
version of the whole, up to a certain recursion depth.
center: the center of the shape
depth: the maximal recursion depth
shrink_factor: the factor by which all lengths decrease for a substitution
subshapes: the set of shapes that a shape is substituted with for a recursion step
"""
struct FractalShape{F<:AbstractFloat,I<:Integer,S<:Shape} <: Shape
    center::Vector{F}
    depth::I
    shrink_factor::F
    subshapes::Vector{S}
end

default_intersection(shape::FractalShape) = default_intersection(first(shape.subshapes))

"""
Create a Menger sponge of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function Menger_sponge(
    center::Vector{F},
    R::AbstractFloat,
    depth::I,
)::FractalShape{F,I,Cube{F}} where {F<:AbstractFloat,I<:Integer}
    subcubes = Cube{F}[]
    R_subcube = R / 3
    R_subcube = convert(F, R_subcube)

    for ordinals in product(1:3, 1:3, 1:3)
        ordinals = collect(ordinals)
        m = countmap(ordinals)
        if 2 ∈ keys(m) && countmap(ordinals)[2] > 1
            continue
        end
        center_subcube = @. center + (ordinals - 2) * 2 * R_subcube
        center_subcube = convert(Vector{F}, center_subcube)

        subcube = Cube(center_subcube, R_subcube)
        push!(subcubes, subcube)
    end
    return FractalShape{F,I,Cube{F}}(center, depth, 3.0, subcubes)
end

"""
A shape consisting of triangles.
vertices: A n_vertices x 3 matrix of vertices
faces: A n_faces x 3 matrix of faces. Each face is defined 
    by the indices of 3 vertices.
center: The center of the shape
n_vertices: The number of vertices
n_faces: the number of faces
convex: Whether the triangles enclose a convex volume.
"""
struct TriangleShape{F<:AbstractFloat,I<:Integer} <: Shape
    vertices::Matrix{F} # (n_vertices, 3)
    faces::Matrix{I} # (n_faces, 3)
    center::Vector{F}
    n_vertices::I
    n_faces::I
    convex::Bool
end

struct Intersection_with_face{F<:AbstractFloat,I<:Integer} <: Intersection
    t::Vector{F}
    face::Vector{I}
end

default_intersection(::TriangleShape{F,I}) where {F,I} =
    Intersection_with_face{F,I}([Inf], [0])

function Tetrahedron(center::Vector{F}, R::F)::TriangleShape{F} where {F<:AbstractFloat}
    vertices = zeros(F, 4, 3)
    vertices[4, :] .= center + [0.0, 0.0, R]

    ϕ = 2π / 3
    for (i, θ) in enumerate(range(0, 2π, 4)[1:end-1])
        @. vertices[i, :] = center + R * [
            cos(θ) * sin(ϕ)
            sin(θ) * sin(ϕ)
            cos(ϕ)
        ]
    end

    faces = collect(transpose(hcat(collect(combinations(1:4, 3))...)))

    return TriangleShape(vertices, faces, center, 4, 4, true)
end

"""
Create a Sierpinski pyramid of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function Sierpinski_pyramid(
    center::Vector{F},
    R::AbstractFloat,
    depth::I,
)::FractalShape{F,I,TriangleShape{F}} where {F<:AbstractFloat,I<:Integer}
    subtetrahedra = TriangleShape{F}[]
    R_subtetrahedron = convert(F, R / 2)

    tetrahedron_main = Tetrahedron(center, R)
    for i = 1:4
        subtetrahedron =
            Tetrahedron(center + tetrahedron_main.vertices[i, :] / 2, R_subtetrahedron)
        push!(subtetrahedra, subtetrahedron)
    end

    return FractalShape(center, depth, convert(F, 2.0), subtetrahedra)
end


"""
Compute the intersection of a ray with a sphere as the smallest
real solution to a quadratic polynomial, if it exists.
"""
function intersect(
    ray::Ray{F},
    sphere::Sphere,
)::Intersection_plain{F} where {F<:AbstractFloat}
    (; loc, dir) = ray
    (; center, Rsq) = sphere

    diff = loc - center
    a = norm(dir)^2
    b = 2 * dot(dir, diff)
    c = norm(diff)^2 - Rsq
    discr = b^2 - 4 * a * c

    t_int = if discr >= 0
        t_int = (-b - sqrt(discr)) / (2 * a)
    else
        Inf
    end

    return Intersection_plain{F}([t_int])
end

"""
Compute the intersection of a ray with a cube by computing the intersections
with each of the 6 face planes and then checking whether the intersection is within the face.
"""
function intersect(
    ray::Ray{F},
    cube::Cube{F},
)::Intersection_with_dim{F,<:Integer} where {F<:AbstractFloat}
    (; loc, dir) = ray
    (; center, R) = cube

    intersection = default_intersection(cube)

    for dim = 1:3
        bound_small = center[dim] - R
        diff_bound_small = bound_small - loc[dim]
        dir_dim_positive = (dir[dim] > 0) # dir_dim = 0 not taken into account

        if diff_bound_small > 0.0
            if dir_dim_positive
                t_int_candidate = diff_bound_small / dir[dim]
            else
                return intersection
            end
        else
            bound_big = center[dim] + R
            diff_bound_big = bound_big - loc[dim]

            if diff_bound_big > 0.0
                if dir_dim_positive
                    t_int_candidate = diff_bound_big / dir[dim]
                else
                    t_int_candidate = -diff_bound_small / dir[dim]
                end
            else
                if dir_dim_positive
                    return intersection
                else
                    t_int_candidate = diff_bound_big / dir[dim]
                end
            end
        end

        if t_int_candidate < intersection.t[1]
            other_dim = 0
            candidate = true

            while candidate && other_dim < 3
                other_dim += 1

                if other_dim !== dim
                    loc_int_other_dim_1 = loc[other_dim] + t_int_candidate * dir[other_dim]
                    if loc_int_other_dim_1 > center[other_dim] + R
                        candidate = false
                        continue
                    elseif loc_int_other_dim_1 < center[other_dim] - R
                        candidate = false
                        continue
                    end
                end
            end

            if candidate
                intersection = Intersection_with_dim([t_int_candidate], [dim])
            end
        end
    end
    return intersection
end

"""
Compute the intersection of a ray with a fractal shape.
This is done recursively until the recursion depth of the fractal shape.
To compute the intersection of a ray with a subshape, the ray location is transformed.
"""
function intersect(
    ray::Ray{F},
    fractal_shape::FractalShape{F};
    current_depth::Int = 0,
)::Intersection where {F<:AbstractFloat}
    (; subshapes, depth, shrink_factor) = fractal_shape

    intersection = default_intersection(first(subshapes))
    subshape_intersect = nothing

    for subshape in subshapes
        intersection_candidate = intersect(ray, subshape)

        if intersection_candidate.t[1] < intersection.t[1]
            subshape_intersect = subshape
            intersection = intersection_candidate
        end
    end

    if !isnothing(subshape_intersect) && current_depth < depth
        intersection = default_intersection(first(subshapes))
        loc_transformed =
            shrink_factor * (ray.loc + fractal_shape.center - subshape_intersect.center)
        ray_transformed = Ray(loc_transformed, ray.dir)
        intersection_candidate =
            intersect(ray_transformed, fractal_shape; current_depth = current_depth + 1)
        intersection_candidate.t[1] /= shrink_factor

        if intersection_candidate.t < intersection.t
            intersection = intersection_candidate
        end
    end

    return intersection
end

"""
Compute the intersection of a ray with a triangle given by
the triangle vertices. This constitutes a 3 variable
linear system solve.
"""
function intersect(
    ray::Ray,
    triangle_vertices::Vector{Vector{F}},
)::F where {F<:AbstractFloat}
    (; loc, dir) = ray

    v1 = triangle_vertices[1] - triangle_vertices[3]
    v2 = triangle_vertices[2] - triangle_vertices[3]

    M = zeros(3, 3)
    M[:, 1] = v1
    M[:, 2] = v2
    M[:, 3] = -dir

    y = loc - triangle_vertices[3]
    sol = M \ y

    return if all(sol .>= 0.0) && (sol[1] + sol[2] <= 1.0)
        sol[3]
    else
        Inf
    end
end

"""
Compute the intersection of a ray with a triangle shape
as the smallest intersection time over all triangles.
"""
function intersect(
    ray::Ray{F},
    shape::TriangleShape{F},
)::Intersection_with_face{F,<:Integer} where {F<:AbstractFloat}
    (; vertices, faces) = shape

    intersection = default_intersection(shape)

    for i = 1:shape.n_faces
        triangle_vertices = [vertices[j, :] for j in faces[i, :]]
        intersection_candidate =
            Intersection_with_face{F,Int64}([intersect(ray, triangle_vertices)], [i])
        if intersection_candidate.t < intersection.t
            intersection = intersection_candidate
        end
    end

    return intersection
end
