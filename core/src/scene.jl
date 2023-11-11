"""
Object for holding the data of all cameras and shapes in a scene.
"""
struct Scene{F<:AbstractFloat}
    cameras::Dict{Symbol,Camera{F}}
    shapes_cube::Dict{Symbol,Cube{F}}
    shapes_fractal_shape::Dict{Symbol,FractalShape{F}}
    shapes_implicit_surface::Dict{Symbol,ImplicitSurface{F}}
    shapes_sphere::Dict{Symbol,Sphere{F}}
    shapes_triangle_shape::Dict{Symbol,TriangleShape{F}}
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
Get a vector of those fields of the scene object which are dictionaries of shapes.
"""
function get_shape_dicts(scene::Scene)::Vector{Dict}
    shape_dicts = Dict[]
    for fieldname in fieldnames(Rays.Scene)
        shape_dict = getfield(scene, fieldname)
        if typeof(shape_dict).parameters[2] <: Shape
            push!(shape_dicts, shape_dict)
        end
    end
    return shape_dicts
end

"""
Remove all shapes from the scene.
"""
function clear_shapes!(scene::Scene)::Nothing
    for shape_dicts in get_shape_dicts(scene)
        empty!(shape_dicts)
    end
    return nothing
end

function Base.show(io::IO, scene::Scene)::Nothing
    println(io, typeof(scene))
    println(io, "* Cameras:")
    for camera in values(scene.cameras)
        println(io, "\t", camera)
    end
    println(io, "\n* Shapes:")
    for shape_dict in get_shape_dicts(scene)
        for shape in values(shape_dict)
            println(io, "\t", shape)
        end
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
    names_camera = keys(scene.cameras)
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
Add a shape to the scene.
"""
function Base.push!(
    scene::Scene{F},
    shape::Shape{F};
    replace::Bool = false,
)::Nothing where {F}
    name = shape.name
    if !replace && name_exists(scene, name)
        name_new = unique_name(name, scene)
        shape = @set shape.name = name_new
        @debug "name $name already used in scene, changed to $name_new."
        name = name_new
    end

    # Put shape in shape dictionary of corresponding type
    for shape_dict in get_shape_dicts(scene)
        shape_type = typeof(shape_dict).parameters[2]
        if shape isa shape_type
            shape_dict[name] = shape
            break
        end
    end
    return nothing
end
