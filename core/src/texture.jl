abstract type Texture{F<:AbstractFloat} end

struct UniformTexture{F} <: Texture{F}
    color::Vector{F}
end

struct IntegerMappingTexture{F} <: Texture{F}
    mapping::Matrix{F}
    variable::Symbol
end

struct ColorFieldTexture{F} <: Texture{F}
    field::VectorField{F}
end

function color!(
    color::AbstractVector{F},
    intersection::Intersection{F},
    texture::UniformTexture,
)::Nothing
    color .= texture.color
    return nothing
end

function color!(
    color::AbstractVector{F},
    intersection::Intersection{F},
    texture::IntegerMappingTexture,
)::nothing

    (; variable, mapping) = texture

    # getfield would lead to runtime dispatch
    if variable == :dim
        value = intersection.dim[1]
    elseif variable == :face
        value = intersection.face[1]
    else
        error("Invalid integer intersection variable \'$variable\'.")
    end

    color .= mapping[:, data_value]
    return nothing
end
