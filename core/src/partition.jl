function get_bounding_box(scene::Scene{F}, name::Symbol)::BoundingBox{F} where {F}
	shape = scene.shapes[name]
	transform = scene.transforms[name]
	return get_bounding_box(shape; transform)
end

function get_bounding_box(shape::Shape{F}; transform::Union{AffineTransform{F}, Nothing} = nothing)::BoundingBox{F} where {F}
	bounding_box = get_bounding_box_(shape)
	if !isnothing(transform)
		bounding_box = transform_bounding_box(bounding_box, transform)
	end
	return bounding_box
end

function get_bounding_box_(shape::Union{Cube{F}, Sphere{F}, ImplicitSurface{F}})::BoundingBox{F} where {F}
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

function transform_bounding_box(bounding_box::BoundingBox{F}, transform::AffineTransform{F})::BoundingBox{F} where {F}
	(; coordinates_min, coordinates_max) = bounding_box
	(; scaling, rotation, translation) = transform

	coordinates = hcat(coordinates_min, coordinates_max)
	vertices = zeros(F, 8, 3)
	for (n, indices) ∈ enumerate(Iterators.product([1, 2], [1, 2], [1, 2]))
		for (i, index) ∈ enumerate(indices)
			vertices[n, i] = coordinates[i, index]
		end
	end

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

	coordinates_min_new = minimum(vertices, dims = 1)[:]
	coordinates_max_new = maximum(vertices, dims = 1)[:]
	return BoundingBox(coordinates_min_new, coordinates_max_new)
end

function partition_scene!(scene::Scene{F})::Nothing where {F}
	(; shapes, partition) = scene

	bounding_boxes::Dict{Symbol, BoundingBox{F}} = Dict(name => get_bounding_box(scene, name) for (name, shape) ∈ shapes)
	@show bounding_boxes

	return nothing
end
