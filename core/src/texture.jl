abstract type Texture{F<:AbstractFloat} end

"""
A texture where the whole shape has the same color.
"""
struct UniformTexture{F} <: Texture{F}
    color::Vector{F}
end

function Base.show(io::IO, texture::UniformTexture)::Nothing
    (; color) = texture
    print(io, "<UniformTexture; color = $color>")
    return nothing
end

"""
A texture where an integer is given a color.
The shape of the mapping is (3, n_colors).
"""
struct IntegerMappingTexture{F} <: Texture{F}
    mapping::Matrix{F}
end

function Base.show(io::IO, texture::IntegerMappingTexture)::Nothing
    (; mapping) = texture
    n_colors = size(mapping)[2]
    print(io, "<IntegerMappingTexture; $n_colors colors>")
    return nothing
end

"""
A texture where each point in 3D space is assigned a color
with the (in place) vector function field!.
"""
struct ColorFieldTexture{F} <: Texture{F}
    field!::VectorField{F}
end

function Base.show(io::IO, texture::ColorFieldTexture)::Nothing
    (; field!) = texture
    print(io, "<ColorFieldTexture; field! = $field!>")
    return nothing
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
function color!(
    color::AbstractVector{F},
    cache_int::AbstractVector{Int},
    cache_float::AbstractVector{F},
    texture::UniformTexture,
)::Nothing where {F}
    color .= view(texture.color, :)
    return nothing
end

"""
Apply the color associated with the given integer intersection variable.
"""
function color!(
    color::AbstractVector{F},
    cache_int::AbstractVector{Int},
    cache_float::AbstractVector{F},
    texture::IntegerMappingTexture,
)::Nothing where {F}
    (; mapping) = texture
    color .= view(mapping, :, cache_int[1])
    return nothing
end

"""
Apply the color given by the color field at the intersection location.
"""
function color!(
    color::AbstractVector{F},
    cache_int::AbstractVector{Int},
    cache_float::AbstractVector{F},
    texture::ColorFieldTexture,
)::Nothing where {F}
    loc_int .= view(dir, :)
    loc_int .*= t[1]
    loc_int .+= view(loc, :)

    texture.field!(color, loc_int)
    return nothing
end
