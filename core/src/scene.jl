"""
Object for holding the data of all cameras and shapes in a scene.
"""
struct Scene{
    F<:AbstractFloat,
    MF<:AbstractMatrix{F},
    MI<:AbstractMatrix{Int},
    VFS<:TypedSubArray{F,MF,Base.Slice{Base.OneTo{Int64}}},
    VIS<:TypedSubArray{Int,MI,Base.Slice{Base.OneTo{Int64}}},
}
    cameras::Dict{Symbol,Camera{F}}
    # Heterogenous values: slow
    shapes::Dict{Symbol,Shape{F}}
    textures::Dict{Symbol,Texture{F}}
    transforms::Dict{Symbol,AffineTransform{F}}
    # Fully typed wrapped functions: fast
    intersectors::Dict{Symbol,Intersector{VFS,VIS}}
    texturers::Dict{Symbol,Texturer{VFS,VIS}}
    # Intersection data of rays with shapes
    intersections::Intersection{F,MF,MI}
    # Partition of scene for efficient intersection computations
    partition::Partition{F,Symbol}
end

"""
Check whether the given name already is a name of a camera or a shape.
"""
function name_exists(scene::Scene, name::Symbol)::Bool
    for fieldname ∈ fieldnames(Scene)
        value = getfield(scene, fieldname)
        if isa(value, Dict) && name in keys(value)
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

function Base.show(io::IO, scene::Scene{F,MF})::Nothing where {F,MF}
    println(io, "Scene (array type $MF):")

    println(io, "* Partition:")
    println(io, "\t", scene.partition, "\n")

    println(io, "* Cameras:")
    for camera ∈ values(scene.cameras)
        println(io, "\t", camera)
    end

    println(io, "\n* Shapes:")
    for shape ∈ values(scene.shapes)
        println(io, "\t", shape)
    end

    println(io, "\n* Transforms:")
    for (name, transform) ∈ scene.transforms
        println(io, "\t", name, ": ", transform)
    end

    println(io, "\n* Textures:")
    for (name, texture) ∈ scene.textures
        println(io, "\t", name, ": ", texture)
    end
end

"""
Get an empty scene.
"""
function Scene(;
    matrix_prototype::MF = zeros(Float32, 3, 3),
)::Scene where {MF<:AbstractMatrix{F} where {F<:AbstractFloat}}
    ST = Base.Slice{Base.OneTo{Int64}}
    F = eltype(MF)
    MI = typeof(similar(matrix_prototype, Int))
    VFS = TypedSubArray{F,MF,ST}
    VIS = TypedSubArray{Int,MI,ST}
    cameras = Dict{Symbol,Camera{F}}()
    shapes = Dict{Symbol,Shape{F}}()
    textures = Dict{Symbol,Texture{F}}()
    transforms = Dict{Symbol,AffineTransform{F}}()
    intersectors = Dict{Symbol,Intersector{VFS,VIS}}()
    texturers = Dict{Symbol,Texturer{VFS,VIS}}()
    intersections = Intersection(nthreads(); matrix_prototype)
    partition = Partition(Vector{PartitionNode{F,Symbol}}())
    return Scene(
        cameras,
        shapes,
        textures,
        transforms,
        intersectors,
        texturers,
        intersections,
        partition,
    )
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
    transform::AffineTransform{F};
    matrix_prototype::MF = zeros(Float32, 3, 3),
)::Intersector where {F<:AbstractFloat,MF<:AbstractMatrix{F}}
    ST = Base.Slice{Base.OneTo{Int64}}
    VFS = TypedSubArray{F,MF,ST}
    MI = typeof(similar(matrix_prototype, Int))
    VIS = TypedSubArray{Int,MI,ST}
    intersector! =
        (
            t,
            ray_loc,
            ray_dir,
            ray_camera_loc,
            ray_camera_dir,
            cache_int,
            cache_float,
            name_intersected,
        ) -> begin
            inverse_transform!(
                ray_loc,
                ray_dir,
                ray_camera_loc,
                ray_camera_dir,
                view(cache_float, 1:3),
                transform,
            )
            if _intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, shape)
                name_intersected[1] = shape.name
                transform_t!(t, transform)
            end
            return nothing
        end
    return Intersector{VFS,VIS}(intersector!)
end

"""
Create the anonymous texturer function
(Done outside Base.push for shapes to avoid runtime dispatch)
"""
function create_texturer(
    texture::Texture{F},
    transform::RayTransform{F};
    matrix_prototype::MF = zeros(Float32, 3, 3),
)::Texturer where {F<:AbstractFloat,MF<:AbstractMatrix{F}}
    ST = Base.Slice{Base.OneTo{Int64}}
    VFS = TypedSubArray{F,MF,ST}
    VI = typeof(similar(matrix_prototype, Int))
    VIS = TypedSubArray{Int,VI,ST}
    return Texturer{VFS,VIS}(
        (color, cache_int, cache_float, ray_camera_loc, ray_camera_dir, t) -> begin
            if texture isa ColorFieldTexture
                vec_temp = view(cache_float, 1:3)
                inverse_transform!(ray_camera_loc, ray_camera_dir, vec_temp, transform)
            end
            color!(
                color,
                cache_int,
                cache_float,
                ray_camera_loc,
                ray_camera_dir,
                t,
                texture,
            )
            return nothing
        end,
    )
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
    transform = scene.transforms[shape_name]
    scene.texturers[shape_name] = create_texturer(texture, transform)
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
    scene::Scene{F,MF},
    shape::Shape{F};
    replace::Bool = false,
    texture::Union{Texture{F},Nothing} = nothing,
    transform::Union{RayTransform{F},Nothing} = nothing,
)::Shape{F} where {F,MF}
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
        color = similar(scene.intersections.cache_float, (3,))
        color .= one(F)
        texture = UniformTexture(color)
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
    texture = scene.textures[shape_name]
    scene.intersectors[shape_name] = create_intersector(shape, transform)
    scene.texturers[shape_name] = create_texturer(texture, transform)
    return nothing
end
