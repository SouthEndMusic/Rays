struct Scene{F<:AbstractFloat}
    cameras::Vector{Camera{F}}
    shapes::Vector{<:Shape{F}}
end

function Scene(; x::F = 0.0f0)::Scene{F} where {F<:AbstractFloat}
    return Scene(Camera{F}[], Shape{F}[])
end

function Base.push!(scene::Scene{F}, camera::Camera{F})::Nothing where {F}
    push!(scene.cameras, camera)
    return nothing
end

function Base.push!(scene::Scene{F}, shape::Shape{F})::Nothing where {F}
    push!(scene.shapes, shape)
    return nothing
end
