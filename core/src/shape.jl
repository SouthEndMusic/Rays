abstract type Shape{F<:AbstractFloat} end

struct Sphere{F} <: Shape{F}
    name::Symbol
    R::F
    Rsq::F
end

Sphere(R::F) where {F} = Sphere(snake_case_name(Sphere), R, R^2)

function Base.show(io::IO, sphere::Sphere)::Nothing
    (; name) = sphere
    print(io, "<Sphere \'$name\'>")
    return nothing
end

struct Cube{F} <: Shape{F}
    name::Symbol
    R::F
end

Cube(R::F) where {F} = Cube(snake_case_name(Cube), R)

function Base.show(io::IO, cube::Cube)::Nothing
    (; name) = cube
    print(io, "<Cube \'$name\'>")
    return nothing
end

"""
A shape where each subshape is substituted by a shrinked
version of the whole, up to a certain recursion depth.
depth: the maximal recursion depth
shrink_factor: the factor by which all lengths decrease for a substitution
subshapes: the set of shapes that a shape is substituted with for a recursion step
"""
struct FractalShape{F,S<:Shape{F},T<:RayTransform{F}} <: Shape{F}
    name::Symbol
    depth::Int
    shape::S
    subshape_transforms::Vector{T}
end

function Base.show(io::IO, fractal_shape::FractalShape)::Nothing
    (; name, subshape_transforms, shape) = fractal_shape
    print(
        io,
        "<FractalShape \'$name\'; $(length(subshape_transforms)) subshapes of $shape>",
    )
end

"""
construct a Menger sponge of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function menger_sponge(
    R::F,
    depth::Int,
)::FractalShape{F,Cube{F},AffineTransform{F,F,Missing,Vector{F}}} where {F}
    subcube_transforms = AffineTransform{F,F,Missing,Vector{F}}[]
    R_subcube = R / 3

    for ordinals in Iterators.product(1:3, 1:3, 1:3)
        ordinals = collect(ordinals)
        if count(x -> x == 2, ordinals) > 1
            continue
        end
        center_subcube = (ordinals .- 2) * 2 * R_subcube
        center_subcube = convert(Vector{F}, center_subcube)

        push!(
            subcube_transforms,
            AffineTransform(F(1 / 3), missing, missing, center_subcube),
        )
    end
    return FractalShape(:menger_sponge, depth, Cube(R), subcube_transforms)
end

"""
A shape consisting of triangles.
vertices: A n_vertices x 3 matrix of vertices
faces: A n_faces x 3 matrix of faces. Each face is defined 
	by the indices of 3 vertices.
n_vertices: The number of vertices
n_faces: the number of faces
convex: Whether the triangles enclose a convex volume.
"""
struct TriangleShape{F} <: Shape{F}
    name::Symbol
    vertices::Matrix{F} # (n_vertices, 3)
    faces::Matrix{Int} # (n_faces, 3)
    normals::Matrix{F} # (n_faces, 3)
    n_vertices::Int
    n_faces::Int
    convex::Bool
end

"""
Construct a triangle shape where the normals are computed automatically from the vertices.
"""
function TriangleShape(
    vertices::Matrix{F},
    faces::Matrix{Int};
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

    return TriangleShape(name, vertices, faces, normals, n_vertices, n_faces, convex)
end

function Base.show(io::IO, triangle_shape::TriangleShape)::Nothing
    (; name, n_vertices, n_faces) = triangle_shape
    print(io, "<TriangleShape \'$name\'; with $n_vertices vertices and $n_faces faces>")
    return nothing
end

"""
Construct a tetrahedron as a triangle shape with given radius.
"""
function Tetrahedron(R::F)::TriangleShape{F} where {F}
    vertices = zeros(F, 4, 3)
    vertices[4, :] .= [0.0, 0.0, R]

    ϕ = 2π / 3
    for (i, θ) in enumerate(range(0, 2π, 4)[1:end-1])
        @. vertices[i, :] = R * [
            cos(θ) * sin(ϕ)
            sin(θ) * sin(ϕ)
            cos(ϕ)
        ]
    end

    faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4]

    return TriangleShape(vertices, faces; convex = true, name = :tetrahedron)
end

"""
Construct a Sierpinski pyramid of given location, size and recursion depth,
with the subshapes array automatically generated.
"""
function sierpinski_pyramid(R::F, depth::Int)::FractalShape{F,TriangleShape{F}} where {F}
    subtetrahedron_transforms = AffineTransform{F,F,Missing,Vector{F}}[]
    tetrahedron = Tetrahedron(R)
    for i ∈ 1:4
        center_subtetrahedron = tetrahedron.vertices[i, :] / 2
        subtetrahedron_transform =
            AffineTransform(F(1 / 2), missing, missing, center_subtetrahedron)
        push!(subtetrahedron_transforms, subtetrahedron_transform)
    end

    return FractalShape(:sierpinski_pyramid, depth, tetrahedron, subtetrahedron_transforms)
end

"""
A shape defined as the level 0 set of the function
	f: R^3 → R 
which is a scalar field assumed to be C1 smooth.
Note: the function f must change sign across the level 0 set
to be detected propery.

name: The name of the same
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
    f::ScalarField{F}
    ∇f!::VF
    R_bound::F
    n_divisions::Int
    tol::F
    itermax::Int
end

"""
Construct an implicit surface with optional rootfinding
parameters.
"""
function ImplicitSurface(
    f::Function;
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
    if isnothing(R_bound)
        R_bound = convert(F, 2.0)
    end
    if isnothing(tol)
        tol = convert(F, 1e-3)
    end

    return ImplicitSurface(
        name,
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

"""
A shape defined as the points a distance r(z) from the z-axis,
with closing disks at z_min and z_max.

name: The name of the shape
r: The distance function from the z-axis
dr: The derivative of r. If not provided, 
	a finite difference derivative is used
r_max: The maximum value of r between z_min and z_max
z_min: The minimum z-value (and the location of a closing disk)
z_max: The maximum z-value (and the location of a closing disk)
n_divisions: The amount of steps used between the 2 intersections
	of the bounding cylinder to find a crossing of r
tol: The tolerance of approximating a zero of r:
	|r - <distance to z-axis>| < tol 
itermax: The maximum amount of Newton iterations used to find 
	a crossing of r along a ray within the specified tolerance
"""
struct RevolutionSurface{F,SF<:Union{ScalarFunc{F},Nothing}} <: Shape{F}
    name::Symbol
    r::ScalarFunc{F}
    dr::SF
    r_max::F
    z_min::F
    z_max::F
    n_divisions::Int
    tol::F
    itermax::Int
end

"""
Construct a revolution surface with optional rootfinding
parameters.
"""
function RevolutionSurface(
    r::Function,
    r_max::F,
    z_min::F,
    z_max::F;
    dr::Union{Function,Nothing} = nothing,
    name::Union{Symbol,Nothing} = nothing,
    itermax::Int = 10,
    n_divisions::Int = 3,
    tol::Union{F,Nothing} = nothing,
)::RevolutionSurface{F} where {F}
    if isnothing(name)
        name = snake_case_name(RevolutionSurface)
    end
    if isnothing(tol)
        tol = convert(F, 1e-3)
    end
    return RevolutionSurface(
        name,
        ScalarFunc{F}(r),
        isnothing(dr) ? nothing : ScalarFunc{F}(dr),
        r_max,
        z_min,
        z_max,
        n_divisions,
        tol,
        itermax,
    )
end

function Base.show(io::IO, implicit_surface::RevolutionSurface)::Nothing
    (; name, r, dr) = implicit_surface
    dr_descr = isnothing(dr) ? "finite difference derivative" : "derivative \'$(dr.obj.x)\'"
    print(io, "<RevolutionSurface \'$name\'; function \'$(r.obj.x)\' and $dr_descr>")
    return nothing
end
