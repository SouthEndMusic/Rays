const Intersector = FunctionWrapper{Nothing,Tuple{Intersection{F}}} where {F}
const Texturer = FunctionWrapper{Nothing,Tuple{Intersection{F}}} where {F}

"""
Object for holding the data of all cameras and shapes in a scene.
"""
struct Scene{F<:AbstractFloat}
    cameras::Dict{Symbol,Camera{F}}
    # Heterogenous values: slow
    shapes::Dict{Symbol,Shape{F}}
    textures::Dict{Symbol,Texture{F}}
    transforms::Dict{Symbol,AffineTransform{F}}
    # Fully typed wrapped functions: fast
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
    empty!(scene.transforms)
    empty!(scene.intersectors)
    empty!(scene.textures)
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
    println(io, "\n* Transforms:")
    for (name, transform) in scene.transforms
        println(io, "\t", name, ": ", transform)
    end
    println(io, "\n* Textures:")
    for (name, texture) in scene.textures
        println(io, "\t", name, ": ", texture)
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
If replace = false, a new name with suffix _i is generated for the camera
with i the smallest integer such that the new name is unique in the scene.
NOTE: This creates a new camera instance which is returned by this function.
"""
function Base.push!(
    scene::Scene{F},
    camera::Camera{F};
    replace::Bool = false,
)::Camera{F} where {F}
    name = camera.name
    if !replace && name_exists(scene, name)
        name_new = unique_name(name, scene)
        camera = @set camera.name = name_new
        @debug "name $name already used in scene, changed to $name_new."
        name = name_new
    end
    scene.cameras[name] = camera
    return camera
end

"""
Create the anonymous intersector function
(Done outside Base.push for shapes to avoid runtime dispatch)
"""
function create_intersector(
    shape::Shape{F},
    transform::AffineTransform{F},
)::Intersector{F} where {F}
    return Intersector{F}(
        intersection -> begin
            inverse_transform!(
                intersection.ray,
                intersection.ray_camera,
                transform;
                intersection.vec_temp,
            )
            if _intersect_ray!(intersection, shape)
                intersection.name_intersected[1] = shape.name
                intersection.t[1] = transform_t(intersection, transform)
            end
        end,
    )
end

"""
Create the anonymous texturer function
(Done outside Base.push for shapes to avoid runtime dispatch)
"""
function create_texturer(texture::Texture{F})::Texturer{F} where {F}
    return Texturer{F}(intersection -> color!(intersection, texture))
end

"""
Set the texture of a shape that is already present in the scene.
"""
function set_texture!(
    scene::Scene{F},
    shape_name::Symbol,
    texture::Texture{F},
)::Nothing where {F}
    if !haskey(scene.shapes, shape_name)
        error("Scene contains no shape with name \'$shape_name\'.")
    end
    scene.textures[shape_name] = texture
    scene.texturers[shape_name] = create_texturer(texture)
    return nothing
end

"""
Add a shape to the scene. Features:
- If replace = false, a new name with suffix _i is generated for the shape
  with i the smallest integer such that the new name is unique in the scene.
  NOTE: This creates a new shape instance which is returned by this function.
- If no texture is given, a default uniform white texture is assigned
"""
function Base.push!(
    scene::Scene{F},
    shape::Shape{F};
    replace::Bool = false,
    texture::Union{Texture{F},Nothing} = nothing,
    transform::Union{RayTransform{F},Nothing} = nothing,
)::Shape{F} where {F}
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
    if isnothing(transform)
        transform = identity_transform(F)
    end
    scene.transforms[name] = transform
    # Define intersector for this shape
    scene.intersectors[name] = create_intersector(shape, transform)

    ## Texturer
    if isnothing(texture)
        texture = UniformTexture(ones(F, 3))
    end
    set_texture!(scene, name, texture)

    return shape
end

"""
Set the transform of a shape alrady present in the scene.
"""
function set_transform!(
    scene::Scene{F},
    shape_name::Symbol,
    transform::AffineTransform{F},
)::Nothing where {F}
    # Redefine intersector for this shape
    if !haskey(scene.shapes, shape_name)
        error("Scene contains no shape with name \'$shape_name\'.")
    end
    shape = scene.shapes[shape_name]
    scene.intersectors[shape_name] = create_intersector(shape, transform)
    return nothing
end
