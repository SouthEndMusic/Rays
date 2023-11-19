abstract type Texture{F<:AbstractFloat} end

"""
A texture where the whole shape has the same color.
"""
struct UniformTexture{F} <: Texture{F}
    color::Vector{F}
end

"""
A texture where an integer variable (:dim or :face) is chosen
and each value of this variable is given its own color.
The shape of the mapping is (3, n_colors).
"""
struct IntegerMappingTexture{F} <: Texture{F}
    mapping::Matrix{F}
    variable::Symbol
end

"""
A texture where each point in 3D space is assigned a color
with the (in place) vector function field!.
"""
struct ColorFieldTexture{F} <: Texture{F}
    field!::VectorField{F}
end

"""
Construct a color field texture.
"""
function ColorFieldTexture(field!::Function; F = Float32)::ColorFieldTexture
    return ColorFieldTexture(VectorField{F}(field!))
end

"""
Apply a uniform color.
"""
function color!(intersection::Intersection{F}, texture::UniformTexture)::Nothing where {F}
    (; color) = intersection
    color .= view(texture.color, :)
    return nothing
end

"""
Apply the color associated with the given integer intersection variable.
"""
function color!(
    intersection::Intersection{F},
    texture::IntegerMappingTexture,
)::Nothing where {F}
    (; variable, mapping) = texture
    (; color, face, dim) = intersection

    # getfield would lead to runtime dispatch
    if variable == :dim
        value = dim[1]
    elseif variable == :face
        value = face[1]
    else
        error("Invalid integer intersection variable \'$variable\'.")
    end

    color .= view(mapping, :, value)
    return nothing
end

"""
Apply the color given by the color field at the intersection location.
"""
function color!(
    intersection::Intersection{F},
    texture::ColorFieldTexture,
)::Nothing where {F}
    (; loc_int, ray, t, color) = intersection
    (; loc, dir) = ray

    loc_int .= view(dir, :)
    loc_int .*= t[1]
    loc_int .+= view(loc, :)

    texture.field!(color, loc_int)
    return nothing
end
