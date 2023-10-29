struct Scene{F<:AbstractFloat}
    cameras::Dict{Symbol,Camera{F}}
    shapes::Dict{Symbol,<:Shape{F}}
end

function Base.show(io::IO, scene::Scene)::Nothing
    println(io, typeof(scene))
    println(io, "* Cameras:")
    for camera in values(scene.cameras)
        println(io, "\t", camera)
    end
    println(io, "\n* Shapes:")
    for shape in values(scene.shapes)
        println(io, "\t", shape)
    end
end

"""
Get an empty scene.
"""
function Scene(; x::F = 0.0f0)::Scene{F} where {F<:AbstractFloat}
    return Scene(Dict{Symbol,Camera{F}}(), Dict{Symbol,Shape{F}}())
end

"""
Extend a name with "_i" where i is the smallest positive Integer
such that the new name does not exist yet in the scene.
"""
function unique_name(name_old::Symbol, scene::Scene)::Symbol
    name_base = string(name_old) * "_"
    names_camera = keys(scene.cameras)
    names_shapes = keys(scene.shapes)
    i = 0
    while true
        i += 1
        name = Symbol(name_base * string(i))
        if name ∉ names_camera && name ∉ names_shapes
            return name
        end
    end
end

"""
Add a camera to the scene.
"""
function Base.push!(scene::Scene{F}, camera::Camera{F};)::Nothing where {F}
    name = camera.name[1]
    if name ∈ keys(scene.cameras) || name ∈ keys(scene.shapes)
        name_new = unique_name(name, scene)
        @debug "name $name already used in scene, changed to $name_new."
        camera.name[1] = name
        name = name_new
    end
    scene.cameras[name] = camera
    return nothing
end

"""
Add a shape to the scene.
"""
function Base.push!(scene::Scene{F}, shape::Shape{F};)::Nothing where {F}
    name = shape.name[1]
    if name ∈ keys(scene.cameras) || name ∈ keys(scene.shapes)
        name_new = unique_name(name, scene)
        @debug "name $name already used in scene, changed to $name_new."
        shape.name[1] = name
        name = name_new
    end
    scene.shapes[name] = shape
    return nothing
end
