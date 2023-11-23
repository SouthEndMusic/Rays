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

function forward_transform!(
    ray_dst::Ray{F},
    ray_src::Ray{F},
    transform::AffineTransform{F},
)::Nothing where {F}
    (; scaling, rotation, translation) = transform
    copyto!(ray_dst.loc, ray_src.loc)
    copyto!(ray_dst.dir, ray_src.dir)

    if !isnothing(rotation)
        ray_dst.loc .= rotation * ray_dst.loc
        ray_dst.dir .= rotation * ray_dst.dir
    end
    if !isnothing(scaling)
        ray_dst.loc .*= scaling
    end
    if !isnothing(translation)
        ray_dst.loc .+= view(translation, :)
    end
    return nothing
end

function inverse_transform!(
    ray_dst::Ray{F},
    ray_src::Ray{F},
    transform::AffineTransform{F},
)::Nothing where {F}
    (; scaling, rotation_inverse, translation) = transform
    copyto!(ray_dst.loc, ray_src.loc)
    copyto!(ray_dst.dir, ray_src.dir)

    if !isnothing(translation)
        ray_dst.loc .-= view(translation, :)
    end
    if !isnothing(scaling)
        ray_dst.loc ./= scaling
    end
    if !isnothing(rotation_inverse)
        ray_dst.loc .= rotation_inverse * ray_dst.loc
        ray_dst.dir .= rotation_inverse * ray_dst.dir
    end
    return nothing
end
