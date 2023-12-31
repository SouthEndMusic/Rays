"""
Get the bounding box of a shape in the scene by its name.
"""
function get_bounding_box(scene::Scene{F}, name::Symbol)::BoundingBox{F} where {F}
    shape = scene.shapes[name]
    transform = scene.transforms[name]
    return get_bounding_box(shape; transform)
end

"""
Get the axes-aligned bounding box of a shape after 
its transform.
"""
function get_bounding_box(
    shape::Shape{F};
    transform::Union{AffineTransform{F},Nothing} = nothing,
)::BoundingBox{F} where {F}
    bounding_box = get_bounding_box_(shape)
    if !isnothing(transform)
        bounding_box = transform_bounding_box(bounding_box, transform)
    end
    return bounding_box
end

function get_bounding_box_(
    shape::Union{Cube{F},Sphere{F},ImplicitSurface{F}},
)::BoundingBox{F} where {F}
    R = shape isa ImplicitSurface ? shape.R_bound : shape.R
    coordinates_min = fill(-R, 3)
    coordinates_max = fill(R, 3)
    return BoundingBox(coordinates_min, coordinates_max)
end

function get_bounding_box_(shape::FractalShape{F})::BoundingBox{F} where {F}
    return get_bounding_box_(shape.shape)
end

function get_bounding_box_(shape::RevolutionSurface{F})::BoundingBox{F} where {F}
    (; z_min, z_max, r_max) = shape
    coordinates_min = F[-r_max, -r_max, z_min]
    coordinates_max = F[r_max, r_max, z_max]
    return BoundingBox(coordinates_min, coordinates_max)
end

function get_bounding_box_(shape::TriangleShape{F})::BoundingBox{F} where {F}
    (; vertices) = shape
    coordinates_min = minimum(vertices, dims = 1)[:]
    coordinates_max = maximum(vertices, dims = 1)[:]
    return BoundingBox(coordinates_min, coordinates_max)
end

"""
Get the axes-aligned bounding box around the given bounding 
box with the given transform applied to it.
"""
function transform_bounding_box(
    bounding_box::BoundingBox{F},
    transform::AffineTransform{F},
)::BoundingBox{F} where {F}
    (; coordinates_min, coordinates_max) = bounding_box
    (; scaling, rotation, translation) = transform

    # Get the coordinates of the vertices of the bounding box
    coordinates = hcat(coordinates_min, coordinates_max)
    vertices = zeros(F, 8, 3)
    for (n, indices) ∈ enumerate(Iterators.product([1, 2], [1, 2], [1, 2]))
        for (i, index) ∈ enumerate(indices)
            vertices[n, i] = coordinates[i, index]
        end
    end

    # Apply the transform to the vertices of the bounding box
    if !ismissing(scaling)
        vertices *= scaling
    end

    if !ismissing(rotation)
        for i ∈ 1:8
            vertices[i, :] = rotation * vertices[i, :]
        end
    end

    if !ismissing(translation)
        for i ∈ 1:8
            vertices[i, :] += translation
        end
    end

    # Get the new bounding box from the transformed vertices
    coordinates_min_new = minimum(vertices, dims = 1)[:]
    coordinates_max_new = maximum(vertices, dims = 1)[:]
    return BoundingBox(coordinates_min_new, coordinates_max_new)
end

"""
Get a partition of the scene in the scene.partition field.
The algorithm starts with a bounding box around all shapes in the scene,
and then builds a tree structure of nested bounding boxes by 
bisecting bounding boxes into 2 based on heuristics for the best splitting plane.
"""
function partition_scene!(
    scene::Scene{F};
    max_shapes_per_node::Int = 3,
    max_depth::Int = 5,
)::Nothing where {F}
    (; shapes, partition) = scene

    # Get the bounding boxes of all shapes in the scene
    bounding_boxes::Dict{Symbol,BoundingBox{F}} =
        Dict(name => get_bounding_box(scene, name) for (name, shape) ∈ shapes)

    # Get the centers of the bounding boxes of all shapes in the scene
    bounding_box_centers::Dict{Symbol,Vector{F}} =
        Dict(name => center(bounding_boxes[name]) for name ∈ keys(shapes))

    # Compute the bounding box around all shapes in the scene
    coordinates_min = [
        minimum(
            bounding_box.coordinates_min[dim] for bounding_box ∈ values(bounding_boxes)
        ) for dim ∈ 1:3
    ]
    coordinates_max = [
        maximum(
            bounding_box.coordinates_max[dim] for bounding_box ∈ values(bounding_boxes)
        ) for dim ∈ 1:3
    ]
    bounding_box_scene = BoundingBox(coordinates_min, coordinates_max)

    # Remove previous partition if it exists
    empty!(partition)

    # Add the first node to the partition
    node_first =
        PartitionNode(bounding_box_scene, Int[], collect(keys(bounding_boxes)), 0, 0)
    push!(partition, node_first)

    # Stack of nodes to process
    nodes_to_process = [1]

    while !isempty(nodes_to_process)
        node_index = pop!(nodes_to_process)
        node = partition[node_index]

        # Compute the splitting dimension as the dimension along which there
        # is the highest variance in bounding box centers
        variances = var([bounding_box_centers[name] for name ∈ node.shape_names])
        dim_split = argmax(variances)

        node = @set node.dim_split = dim_split
        partition[node_index] = node

        # Compute the position of the splitting plane as the median of the bounding box centers
        # of the shapes in this node in the splitting dimension
        coordinate_split =
            median([bounding_box_centers[name][dim_split] for name ∈ node.shape_names])

        # Compute to which of the child nodes each shape belongs
        names_lower = Vector{Symbol}()
        names_higher = Vector{Symbol}()
        for name ∈ node.shape_names
            if bounding_boxes[name].coordinates_min[dim_split] < coordinate_split
                push!(names_lower, name)
            end
            if bounding_boxes[name].coordinates_max[dim_split] > coordinate_split
                push!(names_higher, name)
            end
        end

        depth_new = node.depth + 1

        # The child node with bounding box on the lower side along the splitting dimension
        node_lower = PartitionNode(
            BoundingBox(
                copy(node.bounding_box.coordinates_min),
                copy(node.bounding_box.coordinates_max),
            ),
            Int[],
            names_lower,
            depth_new,
            0,
        )
        node_lower.bounding_box.coordinates_max[dim_split] = coordinate_split
        push!(partition, node_lower)
        node_lower_index = length(partition)
        push!(node.children, node_lower_index)

        # If there are sufficient shapes in the child node and the maximum depth has not been exceeded,
        # add the child node to the stack
        if length(names_lower) > max_shapes_per_node && depth_new < max_depth
            push!(nodes_to_process, node_lower_index)
        end

        # The child node with bounding box on the higher side along the splitting dimension
        node_higher = PartitionNode(
            BoundingBox(
                copy(node.bounding_box.coordinates_min),
                copy(node.bounding_box.coordinates_max),
            ),
            Int[],
            names_higher,
            depth_new,
            0,
        )
        node_higher.bounding_box.coordinates_min[dim_split] = coordinate_split
        push!(partition, node_higher)
        node_higher_index = length(partition)
        push!(node.children, node_higher_index)

        # If there are sufficient shapes in the child node and the maximum depth has not been exceeded,
        # add the child node to the stack
        if length(names_higher) > max_shapes_per_node && depth_new < max_depth
            push!(nodes_to_process, node_higher_index)
        end
    end
    return nothing
end
