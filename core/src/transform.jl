abstract type RayTransform{F<:AbstractFloat} end

struct AffineTransform{
    F<:Union{AbstractFloat,Nothing},
    S<:Union{F,Nothing},
    R<:Union{Matrix{F},Nothing},
    T<:Union{Vector{F},Nothing},
} <: RayTransform{F}
    scaling::S
    rotation::R
    rotation_inverse::R
    translation::T
end

function identity_transform(F)::AffineTransform{F,Nothing,Nothing,Nothing}
    return AffineTransform{F,Nothing,Nothing,Nothing}(nothing, nothing, nothing, nothing)
end

function translation(vector::Vector{F})::AffineTransform{F} where {F}
    return AffineTransform(nothing, nothing, nothing, vector)
end

function inverse_transform!(
    ray::Ray{F},
    ray_camera::Ray{F},
    transform::AffineTransform{F},
)::Nothing where {F}
    (; scaling, rotation_inverse, translation) = transform
    ray.loc .= view(ray_camera.loc, :)
    ray.dir .= view(ray_camera.dir, :)

    if !isnothing(translation)
        ray.loc .-= view(translation, :)
    end
    if !isnothing(scaling)
        ray.loc ./= scaling
    end
    if !isnothing(rotation_inverse)
        ray.loc .= rotation_inverse * ray.loc
        ray.dir .= rotation_inverse * ray.dir
    end
    return nothing
end
