abstract type RayTransform{F<:AbstractFloat} end

"""
A ray transform type that can contain a scaling, a rotation
and a translation.
"""
struct AffineTransform{
    F<:AbstractFloat,
    S<:Union{F,Missing},
    R<:Union{Matrix{F},Missing},
    T<:Union{Vector{F},Missing},
} <: RayTransform{F}
    scaling::S
    rotation::R
    rotation_inverse::R
    translation::T
end

function Base.show(io::IO, transform::AffineTransform{F})::Nothing where {F}
    (; scaling, rotation, translation) = transform
    if transform isa AffineTransform{F,Missing,Missing,Missing}
        print(io, "<AffineTransform; identity>")
    else
        s = "<AffineTransform;"
        print(io, "<AffineTransform; with ")
        if !ismissing(scaling)
            s *= "scaling, "
        end
        if !ismissing(rotation)
            s *= "rotation, "
        end
        if !ismissing(translation)
            s *= "translation, "
        end
        s = s[1:end-2] * ">"
        print(io, s)
    end
    return nothing
end

"""
Get the identity (affine) transform of given type.
"""
function identity_transform(F)::AffineTransform{F,Missing,Missing,Missing}
    return AffineTransform{F,Missing,Missing,Missing}(missing, missing, missing, missing)
end

"""
Get an affine transform which only consists of a translation.
"""
function translation(vector::Vector{F})::AffineTransform{F} where {F}
    return AffineTransform(missing, missing, missing, vector)
end

"""
Get an affine transform which only consists of a rotation around the given
axis by the given amount.
Note: The axis must be a unit vector.
"""
function rotation(axis::Vector{F}, θ::F)::AffineTransform{F} where {F}
    R = zeros(F, 3, 3)
    W = [0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0]
    for i ∈ 1:3
        R[i, i] += 1
    end
    R += sin(θ) * W
    R += (2 * sin(θ / 2)^2) * W^2
    return AffineTransform(missing, R, inv(R), missing)
end

"""
Define the * operator for composing (combining) 
2 affine transformations.
"""
function Base.:*(
    transform_second::AffineTransform{F},
    transform_first::AffineTransform{F},
)::AffineTransform{F} where {F}
    scaling_first = coalesce(transform_first.scaling, one(F))
    scaling_second = coalesce(transform_second.scaling, one(F))
    rotation_first = coalesce(transform_first.rotation, identity_matrix(F))
    rotation_second = coalesce(transform_second.rotation, identity_matrix(F))
    translation_first = coalesce(transform_first.translation, zeros(F, 3))
    translation_second = coalesce(transform_second.translation, zeros(F, 3))

    scaling = scaling_first * scaling_second

    rotation = scaling * rotation_second * rotation_first
    if all(rotation ≈ identity_matrix(F))
        rotation = missing
        rotation_inverse = missing
    else
        rotation_inverse = inv(rotation)
    end

    translation = scaling_second * rotation_second * translation_first + translation_second
    if all(translation .≈ zero(F))
        translation = missing
    end

    if scaling ≈ one(F)
        scaling = missing
    end
    return AffineTransform(scaling, rotation, rotation_inverse, translation)
end

"""
Apply an affine transformation on the source ray and store the result
in the destination ray.
"""
function forward_transform!(
    ray_dst::Ray{F},
    ray_src::Ray{F},
    transform::AffineTransform{F};
    vec_temp::Union{Vector{F},Nothing} = nothing,
)::Nothing where {F}
    (; scaling, rotation, translation) = transform
    if ray_src !== ray_dst
        copyto!(ray_dst.loc, ray_src.loc)
        copyto!(ray_dst.dir, ray_src.dir)
    end

    if !ismissing(rotation)
        mul!(vec_temp, rotation, ray_dst.loc)
        copyto!(ray_dst.loc, temp_vector)
        mul!(vec_temp, rotation, ray_dst.dir)
        copyto!(ray_dst.dir, temp_vector)
    end
    if !ismissing(scaling)
        for i ∈ 1:3
            ray_dst.loc[i] *= scaling
        end
    end
    if !ismissing(translation)
        ray_dst.loc .+= view(translation, :)
    end
    return nothing
end

"""
Apply the inverse of an affine transform on the source ray and store
the results in the destination ray.
"""
function inverse_transform!(
    ray_dst::Ray{F},
    ray_src::Ray{F},
    transform::AffineTransform{F};
    vec_temp::Union{Vector{F},Nothing} = nothing,
)::Nothing where {F}
    (; scaling, rotation_inverse, translation) = transform
    if ray_src !== ray_dst
        copyto!(ray_dst.loc, ray_src.loc)
        copyto!(ray_dst.dir, ray_src.dir)
    end

    if !ismissing(translation)
        ray_dst.loc .-= view(translation, :)
    end
    if !ismissing(scaling)
        for i ∈ 1:3
            ray_dst.loc[i] /= scaling
        end
    end
    if !ismissing(rotation_inverse)
        mul!(vec_temp, rotation_inverse, ray_dst.loc)
        copyto!(ray_dst.loc, vec_temp)
        mul!(vec_temp, rotation_inverse, ray_dst.dir)
        copyto!(ray_dst.dir, vec_temp)
    end
    return nothing
end