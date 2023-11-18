const Intersector = FunctionWrapper{Nothing,Tuple{Intersection{F}}} where {F}
const Texturer = FunctionWrapper{Nothing,Tuple{Intersection{F}}} where {F}

"""
Object for holding the data of all cameras and shapes in a scene.
"""
struct Scene{F<:AbstractFloat}
    cameras::Dict{Symbol,Camera{F}}
    shapes::Dict{Symbol,Shape{F}}
    intersectors::Dict{Symbol,Intersector{F}}
    texturers::Dict{Symbol,Texturer{F}}
end

"""
Check whether the given name already is a name of a camera or a shape.
"""
function name_exists(scene::Scene, name::Symbol)::Bool
    for fieldname in fieldnames(Scene)
        if name in keys(getfield(scene, fieldname))
            return true
        end
    end
    return false
end

"""
Remove all shapes from the scene.
"""
function clear_shapes!(scene::Scene)::Nothing
    empty!(scene.shapes)
    empty!(scene.intersectors)
    empty!(scene.texturers)
    return nothing
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
    return Scene([type{F}() for type in fieldtypes(Scene)]...)
end

"""
Extend a name with "_i" where i is the smallest positive Integer
such that the new name does not exist yet in the scene.
"""
function unique_name(name_old::Symbol, scene::Scene)::Symbol
    name_base = string(name_old) * "_"
    i = 0
    while true
        i += 1
        name = Symbol(name_base * string(i))
        if !name_exists(scene, name)
            return name
        end
    end
end

"""
Add a camera to the scene.
"""
function Base.push!(
    scene::Scene{F},
    camera::Camera{F};
    replace::Bool = false,
)::Nothing where {F}
    name = camera.name
    if !replace && name_exists(scene, name)
        name_new = unique_name(name, scene)
        camera = @set camera.name = name_new
        @debug "name $name already used in scene, changed to $name_new."
        name = name_new
    end
    scene.cameras[name] = camera
    return nothing
end

"""
Create the anonymous intersector function
(Done outside Base.push for shapes to avoid runtime dispatch)
"""
function create_intersector(shape::Shape)::Function
    return intersection ->
        _intersect_ray!(intersection, shape) &&
            (intersection.name_intersected[1] = shape.name)
end

"""
Create the anonymous texturer function
(Done outside Base.push for shapes to avoid runtime dispatch)
"""
function create_texturer(texture::Texture)::Function
    return intersection -> color!(intersection, texture)
end

function set_texture!(
    scene::Scene{F},
    shape_name::Symbol,
    texture::Texture{F},
)::Nothing where {F}
    if !haskey(scene.shapes, shape_name)
        error("Scene contains no shape with name \'$shape_name\'.")
    end
    scene.texturers[shape_name] = Texturer{F}(create_texturer(texture))
    return nothing
end

"""
Add a shape to the scene.
"""
function Base.push!(
    scene::Scene{F},
    shape::Shape{F};
    replace::Bool = false,
    texture::Union{Texture{F},Nothing} = nothing,
)::Nothing where {F}
    ## Name
    name = shape.name
    if !replace && name_exists(scene, name)
        name_new = unique_name(name, scene)
        shape = @set shape.name = name_new
        @debug "name $name already used in scene, changed to $name_new."
        name = name_new
    end
    scene.shapes[name] = shape

    ## Intersector
    # Define intersector for this shape
    scene.intersectors[name] = Intersector{F}(create_intersector(shape))
    # Call the intersector once
    scene.intersectors[name](Intersection(; F))

    ## Texturer
    if isnothing(texture)
        texture = UniformTexture(ones(F, 3))
    end
    set_texture!(scene, name, texture)

    return nothing
end
