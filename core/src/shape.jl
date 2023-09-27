abstract type Shape{F<:AbstractFloat} end

# Intersection types only contain one value per field, but they are vectors to make them mutable.

"""
Get default values for the metadata of the intersection of
a certain shape for when there is no intersection.
"""
struct Sphere{F} <: Shape{F}
    center::Vector{F}
    R::F
    Rsq::F
end

Sphere(center::Vector{F}, R::F) where {F} = Sphere(center, R, R^2)

struct Cube{F} <: Shape{F}
    center::Vector{F}
    R::F
end

"""
A shape where each subshape is substituted by a shrinked
version of the whole, up to a certain recursion depth.
center: the center of the shape
depth: the maximal recursion depth
shrink_factor: the factor by which all lengths decrease for a substitution
subshapes: the set of shapes that a shape is substituted with for a recursion step
"""
struct FractalShape{F,S<:Shape} <: Shape{F}
    center::Vector{F}
    depth::Int
    shrink_factor::F
    subshapes::Vector{S}
end

"""
Create a Menger sponge of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function Menger_sponge(
    center::Vector{F},
    R::AbstractFloat,
    depth::Int,
)::FractalShape{F,Cube{F}} where {F}
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
    return FractalShape{F,Cube{F}}(center, depth, 3.0, subcubes)
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
struct TriangleShape{F} <: Shape{F}
    vertices::Matrix{F} # (n_vertices, 3)
    faces::Matrix{Int} # (n_faces, 3)
    normals::Matrix{F} # (n_faces, 3)
    center::Vector{F}
    n_vertices::Int
    n_faces::Int
    convex::Bool
end

"""
Construct a triangle shape where the normals are computed automatically from the vertices.
"""
function TriangleShape(vertices::Matrix{F}, faces::Matrix{Int}, center::Vector{F}; convex::Bool = false)::TriangleShape{F} where {F}
    n_vertices = size(vertices)[1]
    n_faces = size(faces)[1]
    normals = zeros(F, n_faces, 3)

    @threads for i = 1:n_faces
        u = vertices[faces[i,1],:] - vertices[faces[i,3],:]
        v = vertices[faces[i,2],:] - vertices[faces[i,3],:]
        n = cross(u,v)
        normalize!(n)
        normals[i,:] = n
    end

    return TriangleShape(vertices, faces, normals, center, n_vertices, n_faces, convex)
end

"""
Construct a tetrahedron as a triangle shape with given center and radius.
"""
function Tetrahedron(center::Vector{F}, R::F)::TriangleShape{F} where {F}
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

    return TriangleShape(vertices, faces, center; convex = true)
end

"""
Create a Sierpinski pyramid of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function Sierpinski_pyramid(
    center::Vector{F},
    R::AbstractFloat,
    depth::Int,
)::FractalShape{F,TriangleShape{F}} where {F}
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
