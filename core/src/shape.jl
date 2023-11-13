abstract type Shape{F<:AbstractFloat} end

const ScalarField = FunctionWrapper{F,Tuple{Vector{F}}} where {F}
const VectorField = FunctionWrapper{Nothing,Tuple{Vector{F},Vector{F}}} where {F}

struct Sphere{F} <: Shape{F}
    name::Symbol
    center::Vector{F}
    R::F
    Rsq::F
end

Sphere(center::Vector{F}, R::F) where {F} = Sphere(snake_case_name(Sphere), center, R, R^2)

function Base.show(io::IO, sphere::Sphere)::Nothing
    (; name) = sphere
    print(io, "<Sphere \'$name\'>")
    return nothing
end

struct Cube{F} <: Shape{F}
    name::Symbol
    center::Vector{F}
    R::F
end

Cube(center::Vector{F}, R::F) where {F} = Cube(snake_case_name(Cube), center, R)

function Base.show(io::IO, cube::Cube)::Nothing
    (; name) = cube
    print(io, "<Cube \'$name\'>")
    return nothing
end

"""
A shape where each subshape is substituted by a shrinked
version of the whole, up to a certain recursion depth.
center: the center of the shape
depth: the maximal recursion depth
shrink_factor: the factor by which all lengths decrease for a substitution
subshapes: the set of shapes that a shape is substituted with for a recursion step
"""
struct FractalShape{F,S<:Shape{F}} <: Shape{F}
    name::Symbol
    center::Vector{F}
    depth::Int
    shrink_factor::F
    subshapes::Vector{S}
end

function Base.show(io::IO, fractal_shape::FractalShape)::Nothing
    (; name, subshapes) = fractal_shape
    print(
        io,
        "<FractalShape \'$name\'; $(length(subshapes)) subshapes of type $(eltype(subshapes))>",
    )
end

"""
construct a Menger sponge of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function menger_sponge(
    center::Vector{F},
    R::AbstractFloat,
    depth::Int,
)::FractalShape{F,Cube{F}} where {F}
    subcubes = Cube{F}[]
    R_subcube = R / 3
    R_subcube = convert(F, R_subcube)
    i = 0

    for ordinals in Iterators.product(1:3, 1:3, 1:3)
        i += 1
        ordinals = collect(ordinals)
        if count(x -> x == 2, ordinals) > 1
            continue
        end
        center_subcube = @. center + (ordinals - 2) * 2 * R_subcube
        center_subcube = convert(Vector{F}, center_subcube)

        subcube = Cube(Symbol("subcube_$i"), center_subcube, R_subcube)
        push!(subcubes, subcube)
    end
    return FractalShape{F,Cube{F}}(:menger_sponge, center, depth, 3.0, subcubes)
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
    name::Symbol
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
function TriangleShape(
    vertices::Matrix{F},
    faces::Matrix{Int},
    center::Vector{F};
    convex::Bool = false,
    name::Union{Symbol,Nothing} = nothing,
)::TriangleShape{F} where {F}
    n_vertices = size(vertices)[1]
    n_faces = size(faces)[1]
    normals = zeros(F, n_faces, 3)

    @threads for i ∈ 1:n_faces
        u = vertices[faces[i, 1], :] - vertices[faces[i, 3], :]
        v = vertices[faces[i, 2], :] - vertices[faces[i, 3], :]
        n = cross(u, v)
        normalize!(n)
        normals[i, :] = n
    end

    if isnothing(name)
        name = snake_case_name(TriangleShape)
    end

    return TriangleShape(
        name,
        vertices,
        faces,
        normals,
        center,
        n_vertices,
        n_faces,
        convex,
    )
end

function Base.show(io::IO, triangle_shape::TriangleShape)::Nothing
    (; name, n_vertices, n_faces) = triangle_shape
    print(io, "<TriangleShape \'$name\'; with $n_vertices vertices and $n_faces faces>")
    return nothing
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

    faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4]

    return TriangleShape(vertices, faces, center; convex = true, name = :tetrahedron)
end

"""
Construct a Sierpinski pyramid of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function sierpinski_pyramid(
    center::Vector{F},
    R::AbstractFloat,
    depth::Int,
)::FractalShape{F,TriangleShape{F}} where {F}
    subtetrahedra = TriangleShape{F}[]
    R_subtetrahedron = convert(F, R / 2)

    tetrahedron_main = Tetrahedron(center, R)
    for i ∈ 1:4
        subtetrahedron =
            Tetrahedron(center + tetrahedron_main.vertices[i, :] / 2, R_subtetrahedron)
        push!(subtetrahedra, subtetrahedron)
    end

    return FractalShape(:sierpinski_pyramid, center, depth, convert(F, 2.0), subtetrahedra)
end

"""
Create a new instance of the same shape with a different name
"""
function set_name(shape::S, name_new::Symbol)::S where {S<:Shape}
    fields = [name == :name ? name_new : getfield(shape, name) for name in fieldnames(S)]
    return S(fields...)
end

"""
A shape defined as the level 0 set of the function
	f: R^3 → R 
which is a scalar field assumed to be C1 smooth.
Note: the function f must change sign across the level 0 set
to be detected propery.

center: The center of the shape
f: The scalar field function
∇f!: The (in place) gradient function of f. If not provided,
	finite difference gradients are used
R_bound: The radius of the sphere around 'center' in which
	the shape is assumed to be
n_divisions: The amount of steps used between the 2 intersections
	of the bounding sphere to find a sign change in f
tol: The tolerance of approximating a zero of f:
	|f| < tol 
itermax: The maximum amount of Newton iterations used to find 
	a zero of f along a ray within the specified tolerance
"""
struct ImplicitSurface{F,VF<:Union{VectorField{F},Nothing}} <: Shape{F}
    name::Symbol
    center::Vector{F}
    f::ScalarField{F}
    ∇f!::VF
    R_bound::F
    n_divisions::Int
    tol::F
    itermax::Int
end

function Base.convert(
    ::Type{ImplicitSurface{F,VF}},
    implicit_surface::ImplicitSurface,
) where {F,VF}
    return ImplicitSurface{F,VF}(
        [
            getfield(implicit_surface, fieldname) for
            fieldname in fieldnames(ImplicitSurface)
        ]...,
    )
end

"""
Construct an Implicit surface with optional rootfinding
parameters.
"""
function ImplicitSurface(
    f::Function,
    center::Vector{F};
    ∇f!::Union{Function,Nothing} = nothing,
    R_bound::Union{F,Nothing} = nothing,
    name::Union{Symbol,Nothing} = nothing,
    itermax::Int = 10,
    n_divisions::Int = 3,
    tol::Union{F,Nothing} = nothing,
)::ImplicitSurface{F} where {F}
    if isnothing(name)
        name = snake_case_name(ImplicitSurface)
    end
    if isnothing(center)
        center = zeros(F, 3)
    end
    if isnothing(R_bound)
        R_bound = convert(F, 2.0)
    end
    if isnothing(tol)
        tol = convert(F, 1e-3)
    end

    return ImplicitSurface(
        name,
        center,
        ScalarField{F}(f),
        isnothing(∇f!) ? nothing : VectorField{F}(∇f!),
        R_bound,
        n_divisions,
        tol,
        itermax,
    )
end

function Base.show(io::IO, implicit_surface::ImplicitSurface)::Nothing
    (; name, f, ∇f!) = implicit_surface
    gradient_descr =
        isnothing(∇f!) ? "finite difference gradient" : "gradient \'$(∇f!.obj.x)\'"
    print(io, "<ImplicitSurface \'$name\'; function \'$(f.obj.x)\' and $gradient_descr>")
    return nothing
end

struct RevolutionSurface{F} <: Shape{F}
    name::Symbol
    center::Vector{F}
    r::ScalarFunc{F}
    r_max::F
    z_min::F
    z_max::F
end

function RevolutionSurface(
    r::Function,
    r_max::F,
    z_min::F,
    z_max::F,
    center::Vector{F};
    name::Union{Symbol,Nothing} = nothing,
)::RevolutionSurface{F} where {F}
    if isnothing(name)
        name = snake_case_name(RevolutionSurface)
    end
    return RevolutionSurface(name, center, ScalarFunc{F}(r), r_max, z_min, z_max)
end
