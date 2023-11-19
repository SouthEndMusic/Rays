abstract type AffineTransform{F<:AbstractFloat} end

struct IdentityTransform{F} <: AffineTransform{F} end

struct Translation{F<:AbstractFloat} <: AffineTransform{F}
    vector_back::Vector{F}
end

function get_translation(vector::Vector{F})::Translation{F} where {F}
    return Translation(-vector)
end

function affine_transform!(
    ray::Ray{F},
    ray_camera::Ray{F},
    transform::IdentityTransform,
)::Nothing where {F}
    ray.loc .= view(ray_camera.loc, :)
    ray.dir .= view(ray_camera.loc, :)
    return nothing
end

function affine_transform!(
    ray::Ray{F},
    ray_camera::Ray{F},
    transform::Translation{F},
)::Nothing where {F}
    ray.loc .= view(ray_camera.loc, :)
    ray.loc .+= view(transform.vector_back, :)
    ray.dir .= view(ray_camera.dir, :)
    return nothing
end


