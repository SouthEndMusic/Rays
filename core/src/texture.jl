abstract type Texture{F<:AbstractFloat} end

struct UniformTexture{F} <: Texture{F}
    color::Vector{F}
end

struct IntegerMappingTexture{F} <: Texture{F}
    mapping::Matrix{F}
    variable::Symbol
end

struct ColorFieldTexture{F} <: Texture{F}
    field!::VectorField{F}
end

function color!(intersection::Intersection{F}, texture::UniformTexture)::Nothing where {F}
    (; color) = intersection
    color .= view(texture.color, :)
    return nothing
end

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

function color!(
    intersection::Intersection{F},
    texture::ColorFieldTexture,
)::Nothing where {F}
    (; loc_int, ray, t, color) = intersection
    (; loc, dir) = ray

    loc_int .= view(dir, :)
    loc_int .*= t[1]
    loc_int .+= view(loc, :)

    texture.field(color, loc_int)
end
