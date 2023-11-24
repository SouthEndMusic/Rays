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
    if ray_src !== ray_dst
        copyto!(ray_dst.loc, ray_src.loc)
        copyto!(ray_dst.dir, ray_src.dir)
    end

    if !isnothing(rotation)
        ray_dst.loc .= rotation * ray_dst.loc
        ray_dst.dir .= rotation * ray_dst.dir
    end
    if !isnothing(scaling)
        for i = 1:3
            ray_dst.loc[i] *= scaling
        end
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
    if ray_src !== ray_dst
        copyto!(ray_dst.loc, ray_src.loc)
        copyto!(ray_dst.dir, ray_src.dir)
    end

    if !isnothing(translation)
        ray_dst.loc .-= view(translation, :)
    end
    if !isnothing(scaling)
        for i = 1:3
            ray_dst.loc[i] /= scaling
        end
    end
    if !isnothing(rotation_inverse)
        ray_dst.loc .= rotation_inverse * ray_dst.loc
        ray_dst.dir .= rotation_inverse * ray_dst.dir
    end
    return nothing
end
