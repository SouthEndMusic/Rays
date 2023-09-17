abstract type Shape end

struct Sphere <: Shape
    center::Vector{Float64}
    R::Float64
    Rsq::Float64
end

Sphere(center::Vector{Float64}, R::Float64) = Sphere(center, R, R^2)

struct Cube <: Shape
    center::Vector{Float64}
    R::Float64
end

"""
The Menger sponge struct is similar to the Cube struct but has 2 extra fields:
depth: The recursion depth of the Menger sponge fractal
cubes: The subcubes the largest cube consists of, e.g. 1 recursion step
"""
struct Menger_sponge <: Shape
    center::Vector{Float64}
    R::Float64
    depth::Int
    cubes::Vector{Cube}
end

"""
Create a Menger sponge of given location, size and recursion depth
with the 'cubes' array automatically generated.
"""
function Menger_sponge(center::Vector{Float64}, R::Float64, depth::Int)::Menger_sponge
    cubes = Cube[]
    R_cube = R / 3

    for i = 1:3
        for j = 1:3
            for k = 1:3
                m = countmap([i, j, k])
                if 2 ∈ keys(m) && countmap([i, j, k])[2] > 1
                    continue
                end
                center_cube = zeros(3)
                center_cube += center
                center_cube += @. ([i, j, k] - 2) * 2 * R_cube

                cube = Cube(center_cube, R_cube)
                push!(cubes, cube)
            end
        end
    end
    return Menger_sponge(center, R, depth, cubes)
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
Compute the intersection of a ray with a Menger sponge.
This is done recursively until the recursion depth of the Menger sponge.
To compute the intersection of a ray with a sub-cube, the ray location is transformed.
"""
function intersect(
    ray::Ray,
    menger_sponge::Menger_sponge;
    current_depth::Int = 0,
)::Tuple{Float64,NamedTuple}
    t_int = Inf
    dim_int = 0
    cube_intersect = nothing

    for cube in menger_sponge.cubes
        t_int_candidate, int_metadata = intersect(ray, cube)

        if t_int_candidate < t_int
            cube_intersect = cube
            t_int = t_int_candidate
            dim_int = int_metadata.dim_int
        end
    end

    if !isnothing(cube_intersect) && current_depth < menger_sponge.depth
        t_int = Inf
        loc_transformed = 3 * (ray.loc + menger_sponge.center - cube_intersect.center)
        ray_transformed = Ray(loc_transformed, ray.dir)
        t_int_candidate, int_metadata =
            intersect(ray_transformed, menger_sponge; current_depth = current_depth + 1)
        t_int_candidate /= 3.0
        if t_int_candidate < t_int
            dim_int = int_metadata.dim_int
            t_int = t_int_candidate
        end
    end

    return t_int, (; dim_int)
end