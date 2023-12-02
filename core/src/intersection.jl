"""
Object that holds all information of the intersections of rays
with shapes.
"""
struct Intersection{
	F <: AbstractFloat,
	VF <: AbstractVector{F},
	MF <: AbstractMatrix{F},
	MI <: AbstractMatrix{Int},
}
	# A ray as it comes from the camera
	ray_camera::Ray{MF}
	# A ray as affine transformed for a shape
	ray::Ray{MF}
	# For fractalshape intersections
	ray_transformed::Ray{MF}
	# Intersection time
	t::VF
	# Name of intersected shape
	name_intersected::Vector{Symbol}
	# For integer intersection computations
	cache_int::MI
	# For floating point intersection computations
	cache_float::MF
	# color
	color::MF

end

"""
Get default values for the metadata of the intersection of
a certain shape for when there is no intersection.
"""
function Intersection(
	n_intersections::Int;
	size_cache_int::Int = 1,
	size_cache_float::Int = 12,
	vector_prototype::AbstractVector{F} where {F <: AbstractFloat} = zeros(Float32, 3),
)::Intersection
	return Intersection(
		Ray(n_intersections; vector_prototype), # ray_camera
		Ray(n_intersections; vector_prototype), # ray
		Ray(n_intersections; vector_prototype), # ray_transformed
		similar(vector_prototype, n_intersections), # t
		fill(:none, n_intersections), # name_intersected
		similar(vector_prototype, Int, (n_intersections, size_cache_int)), # int_metadata_int
		similar(vector_prototype, (n_intersections, size_cache_float)), # int_metadata_float
		similar(vector_prototype, (n_intersections, 3)), # color
	)
end

"""
Set the intersection to a non-intersected state.
"""
function reset_intersection!(intersection::Intersection, index::Int)::Nothing
	intersection.t[index] = Inf
	intersection.name_intersected[index] = :none
	return nothing
end

"""
Compute the intersections of a ray with a sphere.
Returns (nothing, nothing) if the intersections do not exist.
"""
function intersect_sphere(
	ray::Ray{F},
	Rsq::F,
)::Tuple{Union{F, Nothing}, Union{F, Nothing}} where {F}
	(; loc, dir) = ray
	a = norm(dir)^2
	b = 2 * dot(dir, loc)
	c = norm(loc)^2 - Rsq
	discr = b^2 - 4 * a * c

	if discr >= 0
		denom = 2 * a
		sqrt_discr = sqrt(discr)
		t_int_1 = (-b - sqrt_discr) / denom
		t_int_2 = (-b + sqrt_discr) / denom
		return t_int_1, t_int_2
	else
		return nothing, nothing
	end
end


"""
Compute the intersection of a ray with a sphere as the smallest
real solution to a quadratic polynomial, if it exists.
"""
function _intersect_ray!(intersection::Intersection{F}, sphere::Sphere)::Bool where {F}
	(; ray) = intersection
	(; Rsq) = sphere
	t_int_candidate = intersect_sphere(ray, Rsq)[1]

	closer_intersection_found = false

	if !isnothing(t_int_candidate)
		if t_int_candidate < intersection.t[1]
			closer_intersection_found = true
			intersection.t[1] = t_int_candidate
		end
	end

	return closer_intersection_found
end

"""
Compute the intersection of a ray with a cube by computing the intersections
with each of the 6 face planes and then checking whether the intersection is within the face.
"""
function _intersect_ray!(
	t::AbstractVector{F},
	cache_int::AbstractVector{Int},
	cache_float::AbstractVector{F},
	ray_loc::AbstractVector{F},
	ray_dir::AbstractVector{F},
	cube::Cube{F},
)::Bool where {F}
	closer_intersection_found = false

	for dim ∈ 1:3
		diff_bound_small = -cube.R - ray_loc[dim]
		dir_dim_positive = (ray_dir[dim] > 0) # dir_dim = 0 not taken into account

		if diff_bound_small > 0.0
			if dir_dim_positive
				t_int_candidate = diff_bound_small / ray_dir[dim]
			else
				return closer_intersection_found
			end
		else
			diff_bound_big = cube.R - ray_loc[dim]

			if diff_bound_big > 0.0
				if dir_dim_positive
					t_int_candidate = diff_bound_big / ray_dir[dim]
				else
					t_int_candidate = -diff_bound_small / ray_dir[dim]
				end
			else
				if dir_dim_positive
					return closer_intersection_found
				else
					t_int_candidate = diff_bound_big / ray_dir[dim]
				end
			end
		end

		if t_int_candidate < t[1]
			other_dim = 0
			candidate = true

			while candidate && other_dim < 3
				other_dim += 1

				if other_dim !== dim
					loc_int_other_dim_1 =
						ray_loc[other_dim] + t_int_candidate * ray_dir[other_dim]
					if loc_int_other_dim_1 > cube.R
						candidate = false
						continue
					elseif loc_int_other_dim_1 < -cube.R
						candidate = false
						continue
					end
				end
			end

			if candidate
				closer_intersection_found = true
				t[1] = t_int_candidate
				cache_int[1] = dim
			end
		end
	end
	return closer_intersection_found
end

function _intersect_ray!(
	intersection::Intersection{F},
	fractal_shape::FractalShape{F, S, T};
	current_depth::Int = 1,
)::Bool where {F, S, T}
	(; ray) = intersection
	(; depth, shape, subshape_transforms) = fractal_shape
	closer_intersection_found = false
	reset_intersection!(intersection)

	for subshape_transform in subshape_transforms
		intersection.t[1] /= subshape_transform.scaling
		inverse_transform!(ray, ray, subshape_transform)
		if _intersect_ray!(intersection, shape)
			if current_depth < depth
				if _intersect_ray!(
					intersection,
					fractal_shape,
					current_depth = current_depth + 1,
				)
					closer_intersection_found = true
				end
			else
				closer_intersection_found = true
			end
		end
		forward_transform!(ray, ray, subshape_transform)
		intersection.t[1] *= subshape_transform.scaling
	end
	return closer_intersection_found
end

"""
Compute the intersection of a ray with a triangle given by
the triangle vertices.
"""
function _intersect_ray!(
	intersection::Intersection{F},
	triangle_vertices::AbstractArray{F};
	normal::Union{AbstractVector{F}, Nothing} = nothing,
	t_int_prev::Union{F, Nothing} = nothing,
)::F where {F}
	(; ray, diff, diff2, u, v) = intersection
	(; loc, dir) = intersection.ray

	if isnothing(normal)
		c = view(triangle_vertices, 3, :)
		u .= view(triangle_vertices, 1, :) - c
		v .= view(triangle_vertices, 2, :) - c
		normal = cross(u, v)
		normalize!(normal)
	end

	diff .= 0.0
	diff .+= view(triangle_vertices, 3, :)
	diff .-= loc
	t_int_candidate = dot(diff, normal) / dot(dir, normal)

	if t_int_candidate < 0.0
		return t_int_prev
	end

	if isnothing(t_int_prev)
		t_int_prev = convert(F, Inf)
	end

	if t_int_candidate ≥ t_int_prev
		return t_int_prev
	end

	if !isnothing(normal)
		c = view(triangle_vertices, 3, :)
		u .= 0.0
		v .= 0.0
		u .+= view(triangle_vertices, 1, :)
		v .+= view(triangle_vertices, 2, :)
		u .-= c
		v .-= c
	end

	u_normsq = dot(u, u)
	v_normsq = dot(v, v)
	inner_uv = dot(u, v)
	det = u_normsq * v_normsq - inner_uv^2
	diff2 .= @. t_int_candidate * ray.dir - diff
	inner_u = dot(diff2, u)
	inner_v = dot(diff2, v)
	λ_1 = (inner_u * v_normsq - inner_v * inner_uv) / det

	if !(0 ≤ λ_1 ≤ 1)
		return t_int_prev
	end

	λ_2 = (inner_v * u_normsq - inner_u * inner_uv) / det

	if !(0 ≤ λ_2 ≤ 1)
		return t_int_prev
	end

	if λ_1 + λ_2 ≤ 1.0
		return t_int_candidate
	else
		return t_int_prev
	end
end

"""
Compute the intersection of a ray with a triangle shape
as the smallest intersection time over all triangles.
"""
function _intersect_ray!(
	intersection::Intersection{F},
	shape::TriangleShape{F};
)::Bool where {F}
	(; vertices, faces, normals) = shape

	closer_intersection_found = false

	for i ∈ 1:shape.n_faces
		triangle_vertices = view(vertices, view(faces, i, :), :)
		normal = view(normals, i, :)
		t_int_candidate = _intersect_ray!(
			intersection,
			triangle_vertices;
			normal,
			t_int_prev = intersection.t[1],
		)
		if t_int_candidate < intersection.t[1]
			closer_intersection_found = true
			intersection.t[1] = t_int_candidate
			intersection.face[1] = i
		end
	end

	return closer_intersection_found
end

"""
Compute the gradient of a scalar field with finite differences.
"""
function ∇f_finitediff!(
	grad::Vector{F},
	loc::Vector{F},
	f::ScalarField{F};
	eps::AbstractFloat = 1e-4,
) where {F}
	for i ∈ 1:3
		loc_perturbed = copy(loc)
		loc_perturbed[i] += eps
		grad[i] = (f(loc_perturbed) - f(loc)) / eps
	end
	return nothing
end

"""
Compute the intersection between a ray and an implicit surface.
"""
function _intersect_ray!(
	intersection::Intersection{F},
	shape::ImplicitSurface{F},
)::Bool where {F}
	(; loc_int, ray, grad) = intersection
	(; loc, dir) = ray
	(; f, ∇f!, itermax, tol, R_bound, n_divisions) = shape

	# Compute intersections of the ray with the bounding sphere
	bound_lower, bound_upper = intersect_sphere(ray, R_bound^2)

	# If there are no intersections with the bounding sphere,
	# there are certainly no intersections with the implicit surface
	if isnothing(bound_lower)
		return false
	end

	# Compute the value of f at the closer bounding sphere intersection
	loc_int .= view(dir, :)
	loc_int .*= bound_lower
	loc_int .+= view(loc, :)
	fval_lower = f(loc_int)

	# Find a sign change in f between the closer and further
	# bounding sphere intersections starting from the closer
	# and stepping with Δt
	t_0 = zero(F)
	Δt = (bound_upper - bound_lower) / n_divisions
	found_t_0 = false
	for i ∈ 1:n_divisions
		t_0 = bound_lower + i * Δt
		loc_int .= view(dir, :)
		loc_int .*= t_0
		loc_int .+= view(loc, :)
		if fval_lower * f(loc_int) < zero(F)
			found_t_0 = true
			break
		end
	end

	# If no sign change was found, conclude that 
	# there is no intersection of this implicit surface
	# with this ray
	if !found_t_0
		return false
	end

	# Find the zero of f along the ray with the desired
	# tolerance using Newton iterations
	t_n = t_0
	for i ∈ 1:itermax
		loc_int .= view(dir, :)
		loc_int *= t_n
		loc_int += view(loc, :)
		fval = f(loc_int)
		error = abs(fval)
		if error < tol
			if t_n < intersection.t[1]
				intersection.t[1] = t_n
				return true
			else
				return false
			end
		else
			# If no analytical gradient of f is provided,
			# compute a finite difference gradient.
			if isnothing(∇f!)
				∇f_finitediff!(grad, loc_int, f)
			else
				∇f!(grad, loc_int)
			end
			t_n -= fval / dot(grad, dir)
		end
	end
	return false
end

function _intersect_ray!(
	intersection::Intersection{F},
	shape::RevolutionSurface{F},
)::Bool where {F}

	(; loc_int, ray) = intersection
	(; loc, dir) = ray
	(; z_min, z_max, r, r_max, itermax, tol, n_divisions, dr) = shape

	closer_intersection_found = false

	if dir[3] ≈ zero(F)
		return closer_intersection_found
	end

	# Top disk
	t_top = (z_max - loc[3]) / dir[3]
	r_top = sqrt((loc[1] + dir[1] * t_top)^2 + (loc[2] + dir[2] * t_top)^2)
	if r_top <= r(z_max)
		if t_top < intersection.t[1]
			closer_intersection_found = true
			intersection.t[1] = t_top
		end
	end

	# Bottom disk
	t_bottom = (z_min - loc[3]) / dir[3]
	r_bottom = sqrt((loc[1] + dir[1] * t_bottom)^2 + (loc[2] + dir[2] * t_bottom)^2)
	if r_bottom <= r(z_min)
		if t_bottom < intersection.t[1]
			closer_intersection_found = true
			intersection.t[1] = t_bottom
		end
	end

	# Bounding cylinder
	a = dir[1]^2 + dir[2]^2
	b = 2 * (dir[1] * loc[1] + dir[2] * loc[2])
	c = loc[1]^2 + loc[2]^2
	discr = b^2 - 4 * a * (c - r_max^2)
	if discr < zero(F)
		return closer_intersection_found
	end

	denom = 2 * a
	sqrt_discr = sqrt(discr)
	t_int_min = (-b - sqrt_discr) / denom
	t_int_max = (-b + sqrt_discr) / denom

	if dir[3] > 0
		t_min = max(t_int_min, t_bottom)
		t_max = min(t_int_max, t_top)
	else
		t_min = max(t_int_min, t_top)
		t_max = min(t_int_max, t_bottom)
	end


	# Find a sign change from negative to positive in
	# f(t) = r(o_z + d_z*t) - sqrt((o_x+d_x*t)^2 + (o_y+d_y*t)^2) 
	# between the closer and further
	# bounding sphere intersections starting from the closer
	# and stepping with Δt
	t_0 = zero(F)
	Δt = (t_max - t_min) / n_divisions
	found_t_0 = false
	for i ∈ 1:n_divisions
		t_0 = t_min + i * Δt
		fval = r(loc[3] + dir[3] * t_0) - sqrt(a * t_0^2 + b * t_0 + c)
		if fval >= 0
			found_t_0 = true
			break
		end
	end

	if !found_t_0
		return closer_intersection_found
	end

	t_n = t_0
	for i ∈ 1:itermax
		loc_int .= view(dir, :)
		loc_int *= t_n
		loc_int += view(loc, :)
		dist = sqrt(a * t_n^2 + b * t_n + c)
		z_n = loc[3] + dir[3] * t_n
		fval = r(z_n) - dist
		error = abs(fval)
		if error < tol
			if t_n < intersection.t[1]
				intersection.t[1] = t_n
				closer_intersection_found = true
			end
			break
		else
			if isnothing(dr)
				eps = 1e-4
				dr_value = (r(z_n + eps) - r(z_n)) / eps
			else
				dr_value = dr(z_n)
			end
			df = dir[3] * dr_value - (b / 2 + a * t_n) / dist
			t_n -= fval / df
			if !(t_min <= t_n <= t_max)
				break
			end
		end
	end

	return closer_intersection_found
end

"""
Transform the intersection time by the scaling of the given affine transform.
"""
function transform_t!(
	t::AbstractVector{F},
	transform::AffineTransform{F},
)::Nothing where {F}
	if !ismissing(transform.scaling)
		t[1] *= transform.scaling
	end
	return nothing
end
