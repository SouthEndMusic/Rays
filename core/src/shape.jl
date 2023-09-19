abstract type Shape end
abstract type TriangleShape <: Shape end

"""
Get default values for the metadata of the intersection of
a certain shape for when there is no intersection.
"""
default_metadata(shape_type::Shape) = default_metadata(Val(typeof(shape_type)))

struct Sphere <: Shape
    center::Vector{Float64}
    R::Float64
    Rsq::Float64
end

default_metadata(::Val{Sphere}) = (;)
Sphere(center::Vector{Float64}, R::Float64) = Sphere(center, R, R^2)

struct Cube <: Shape
    center::Vector{Float64}
    R::Float64
end

default_metadata(::Val{Cube}) = (; dim_int = 0)

"""
A shape where each subshape is substituted by a shrinked
version of the whole, up to a certain recursion depth.
center: the center of the shape
depth: the maximal recursion depth
shrink_factor: the factor by which all lengths decrease for a substitution
subshapes: the set of shapes that a shape is substituted with for a recursion step
"""
struct FractalShape{T<:Shape} <: Shape
    center::Vector{Float64}
    depth::Int
    shrink_factor::Float64
    subshapes::Vector{T}
end

"""
Create a Menger sponge of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function Menger_sponge(center::Vector{Float64}, R::Float64, depth::Int)::FractalShape
    subcubes = Cube[]
    R_subcube = R / 3

    for i = 1:3
        for j = 1:3
            for k = 1:3
                m = countmap([i, j, k])
                if 2 ∈ keys(m) && countmap([i, j, k])[2] > 1
                    continue
                end
                center_subcube = zeros(3)
                center_subcube += center
                center_subcube += @. ([i, j, k] - 2) * 2 * R_subcube

                subcube = Cube(center_subcube, R_subcube)
                push!(subcubes, subcube)
            end
        end
    end
    return FractalShape(center, depth, 3.0, subcubes)
end

"""
A shape consisting of triangles.
No assumptions are made about the structure of this shape.
vertices: A n_vertices x 3 matrix of vertices
faces: A n_faces x 3 matrix of faces. Each face is defined 
    by the indices of 3 vertices.
center: The center of the shape
n_vertices: The number of vertices
n_faces: the number of faces
"""
struct GeneralTriangleShape <: TriangleShape
    vertices::Matrix{Float64} # (n_vertices, 3)
    faces::Matrix{Int} # (n_faces, 3)
    center::Vector{Float64}
    n_vertices::Int
    n_faces::Int
end

default_metadata(::Val{GeneralTriangleShape}) = (; face_int = 0)

"""
See GeneralTriangleShape. 
This struct assumes that the triangles enclose a convex volume.
"""
struct ConvexTriangleShape <: TriangleShape
    vertices::Matrix{Float64} # (n_vertices, 3)
    faces::Matrix{Int} # (n_faces, 3)
    center::Vector{Float64}
    n_vertices::Int
    n_faces::Int
end

default_metadata(::Val{ConvexTriangleShape}) = (; face_int = 0)

function Tetrahedron(center::Vector{Float64}, R::Float64)::ConvexTriangleShape
    vertices = zeros(4, 3)
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

    return ConvexTriangleShape(vertices, faces, center, 4, 4)
end

"""
Create a Sierpinski pyramid of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function Sierpinski_pyramid(center::Vector{Float64}, R::Float64, depth::Int)::FractalShape
    subtetrahedra = ConvexTriangleShape[]
    R_subtetrahedron = R / 2

    tetrahedron_main = Tetrahedron(center, R)
    for i = 1:4
        subtetrahedron =
            Tetrahedron(center + tetrahedron_main.vertices[i, :] / 2, R_subtetrahedron)
        push!(subtetrahedra, subtetrahedron)
    end

    return FractalShape(center, depth, 2.0, subtetrahedra)
end


"""
Compute the intersection of a ray with a sphere as the smallest
real solution to a quadratic polynomial, if it exists.
"""
function intersect(ray::Ray, sphere::Sphere)::Tuple{Float64,NamedTuple}
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

    return t_int, (;)
end

"""
Compute the intersection of a ray with a cube by computing the intersections
with each of the 6 face planes and then checking whether the intersection is within the face.
"""
function intersect(ray::Ray, cube::Cube)::Tuple{Float64,NamedTuple}
    (; loc, dir) = ray
    (; center, R) = cube

    diff = center - loc

    t_int = Inf
    dim_int = 0

    for r in [-R, R]
        for dim = 1:3
            t_int_candidate = (diff[dim] + r) / dir[dim]
            valid_candidate = true
            for dim_other = 1:3
                if dim == dim_other
                    continue
                end
                loc_intersect = loc[dim_other] + t_int_candidate * dir[dim_other]
                if abs(loc_intersect - center[dim_other]) > R
                    valid_candidate = false
                    break
                end
            end
            if valid_candidate
                if t_int_candidate < t_int
                    dim_int = dim
                    t_int = t_int_candidate
                end
            end
        end
    end
    return t_int, (; dim_int)
end

"""
Compute the intersection of a ray with a fractal shape.
This is done recursively until the recursion depth of the fractal shape.
To compute the intersection of a ray with a subshape, the ray location is transformed.
"""
function intersect(
    ray::Ray,
    fractal_shape::FractalShape;
    current_depth::Int = 0,
)::Tuple{Float64,NamedTuple}
    (; subshapes, depth, shrink_factor) = fractal_shape

    t_int = Inf
    metadata_int = default_metadata(first(subshapes))
    subshape_intersect = nothing

    for subshape in subshapes
        t_int_candidate, metadata_int_candidate = intersect(ray, subshape)

        if t_int_candidate < t_int
            subshape_intersect = subshape
            t_int = t_int_candidate
            metadata_int = metadata_int_candidate
        end
    end

    if !isnothing(subshape_intersect) && current_depth < depth
        t_int = Inf
        loc_transformed =
            shrink_factor * (ray.loc + fractal_shape.center - subshape_intersect.center)
        ray_transformed = Ray(loc_transformed, ray.dir)
        t_int_candidate, metadata_int_candidate =
            intersect(ray_transformed, fractal_shape; current_depth = current_depth + 1)
        t_int_candidate /= shrink_factor
        if t_int_candidate < t_int
            t_int = t_int_candidate
            metadata_int = metadata_int_candidate
        end
    end

    return t_int, metadata_int
end

"""
Compute the intersection of a ray with a triangle given by
the triangle vertices. This constitutes a 3 variable
linear system solve.
"""
function intersect(ray::Ray, triangle_vertices::Vector{Vector{Float64}})::Float64
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
function intersect(ray::Ray, shape::TriangleShape)::Tuple{Float64,NamedTuple}
    (; vertices, faces) = shape

    t_int = Inf
    face_int = 0

    for i = 1:shape.n_faces
        triangle_vertices = [vertices[j, :] for j in faces[i, :]]
        t_int_candidate = intersect(ray, triangle_vertices)
        if t_int_candidate < t_int
            t_int = t_int_candidate
            face_int = i
        end
    end

    return t_int, (; face_int)
end
