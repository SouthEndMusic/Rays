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

    # Get the coordinates of the vertices of the bounding box
    coordinates = hcat(coordinates_min, coordinates_max)
    vertices = zeros(F, 8, 3)
    for (n, indices) ∈ enumerate(Iterators.product([1, 2], [1, 2], [1, 2]))
        for (i, index) ∈ enumerate(indices)
            vertices[n, i] = coordinates[i, index]
        end
    end

    # Apply the transform to the vertices of the bounding box
    forward_transform!(vertices, transform)

    # Get the new bounding box from the transformed vertices
    coordinates_min_new = minimum(vertices, dims = 1)[:]
    coordinates_max_new = maximum(vertices, dims = 1)[:]
    return BoundingBox(coordinates_min_new, coordinates_max_new)
end

"""
Get a partition of the given bounding boxes.
The algorithm starts with a bounding box around all objects,
and then builds a tree structure of nested bounding boxes by 
bisecting bounding boxes into 2 based on heuristics for the best splitting plane.
partition: vector of nodes in the bisection tree
bounding_boxes: dictionary identifier => bounding_box
max_objects_per_node: stop bisecting if a bounding box contains at most this amount
	of objects
max_depth: maximum depth of the bisection tree
"""
function partition!(
    partition::Vector{PartitionNode{F,T}},
    bounding_boxes::Dict{T,BoundingBox{F}};
    max_objects_per_node::Int = 3,
    max_depth::Int = 5,
)::Nothing where {F,T}

    # Get the centers of all bounding boxes
    bounding_box_centers::Dict{T,Vector{F}} = Dict(
        identifier => center(bounding_box) for (identifier, bounding_box) ∈ bounding_boxes
    )

    # Compute the bounding box around all objects
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
    bounding_box_outter = BoundingBox(coordinates_min, coordinates_max)

    # Remove previous partition if it exists
    empty!(partition)

    # Add the first node to the partition
    identifiers = collect(keys(bounding_boxes))
    node_first = PartitionNode(bounding_box_outter, Int[], identifiers, 0, 0)
    push!(partition, node_first)

    # Stack of nodes to process
    nodes_to_process = [1]

    while !isempty(nodes_to_process)
        node_index = pop!(nodes_to_process)
        node = partition[node_index]

        if isempty(node.identifiers)
            continue
        end

        # Compute the splitting dimension as the dimension along which there
        # is the highest variance in bounding box centers
        variances = var([bounding_box_centers[name] for name ∈ node.identifiers])
        dim_split = argmax(variances)

        node = @set node.dim_split = dim_split
        partition[node_index] = node

        # Compute the position of the splitting plane as the median of the bounding box centers
        # of the shapes in this node in the splitting dimension
        coordinate_split =
            median([bounding_box_centers[name][dim_split] for name ∈ node.identifiers])

        # Compute to which of the child nodes each object belongs
        identifiers_lower = Vector{T}()
        identifiers_higher = Vector{T}()
        for identifier ∈ node.identifiers
            if bounding_boxes[identifier].coordinates_min[dim_split] < coordinate_split
                push!(identifiers_lower, identifier)
            end
            if bounding_boxes[identifier].coordinates_max[dim_split] > coordinate_split
                push!(identifiers_higher, identifier)
            end
        end

        depth_new = node.depth + 1

        # The child node with bounding box on the lower side along the splitting dimension
        if length(identifiers_lower) > 0
            node_lower = PartitionNode(
                BoundingBox(
                    copy(node.bounding_box.coordinates_min),
                    copy(node.bounding_box.coordinates_max),
                ),
                Int[],
                identifiers_lower,
                depth_new,
                0,
            )
            node_lower.bounding_box.coordinates_max[dim_split] = coordinate_split
            push!(partition, node_lower)
            node_lower_index = length(partition)
            push!(node.child_indices, node_lower_index)

            # If there are sufficient shapes in the child node and the maximum depth has not been exceeded,
            # add the child node to the stack
            if length(identifiers_lower) > max_objects_per_node && depth_new < max_depth
                push!(nodes_to_process, node_lower_index)
            end
        end

        # The child node with bounding box on the higher side along the splitting dimension
        if length(identifiers_higher) > 0
            node_higher = PartitionNode(
                BoundingBox(
                    copy(node.bounding_box.coordinates_min),
                    copy(node.bounding_box.coordinates_max),
                ),
                Int[],
                identifiers_higher,
                depth_new,
                0,
            )
            node_higher.bounding_box.coordinates_min[dim_split] = coordinate_split
            push!(partition, node_higher)
            node_higher_index = length(partition)
            push!(node.child_indices, node_higher_index)

            # If there are sufficient shapes in the child node and the maximum depth has not been exceeded,
            # add the child node to the stack
            if length(identifiers_higher) > max_objects_per_node && depth_new < max_depth
                push!(nodes_to_process, node_higher_index)
            end
        end
    end
    return nothing
end

function partition!(
    scene::Scene{F};
    max_objects_per_node::Int = 3,
    max_depth::Int = 5,
)::Nothing where {F}
    (; shapes, partition) = scene

    # Get the bounding boxes of all shapes in the scene
    bounding_boxes::Dict{Symbol,BoundingBox{F}} =
        Dict(name => get_bounding_box(scene, name) for (name, shape) ∈ shapes)

    partition!(partition, bounding_boxes; max_objects_per_node, max_depth)
    return nothing
end

function partition!(
    shape::TriangleShape{F};
    max_objects_per_node::Int = 3,
    max_depth::Int = 5,
)::Nothing where {F}
    (; vertices, faces, partition, n_faces) = shape

    # First collect the bouding boxes in a vector,
    # which is thread-safe
    bounding_boxes_vector = Vector{BoundingBox}(undef, n_faces)
    @batch for face_index ∈ 1:n_faces
        triangle_vertices = view(vertices, view(faces, face_index, :), :)
        bounding_box = BoundingBox(
            minimum(triangle_vertices, dims = 1)[:],
            maximum(triangle_vertices, dims = 1)[:],
        )

        bounding_boxes_vector[face_index] = bounding_box
    end

    bounding_boxes =
        Dict(i => bounding_box for (i, bounding_box) ∈ enumerate(bounding_boxes_vector))
    partition!(partition, bounding_boxes; max_objects_per_node, max_depth)
    return nothing
end

function partition!(
    scene::Scene,
    name::Symbol;
    max_objects_per_node::Int = 3,
    max_depth::Int = 5,
)::Nothing
    shape = scene.shapes[name]
    partition!(shape; max_objects_per_node, max_depth)
    return nothing
end
