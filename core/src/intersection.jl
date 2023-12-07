"""
Object that holds all information of the intersections of rays
with shapes.
"""
struct Intersection{F <: AbstractFloat, MF <: AbstractMatrix{F}, MI <: AbstractMatrix{Int}}
	# A ray as it comes from the camera
	ray_camera::Ray{MF}
	# A ray as affine transformed for a shape
	ray::Ray{MF}
	# For fractalshape intersections
	ray_transformed::Ray{MF}
	# Intersection time
	t::MF
	# Name of intersected shape
	name_intersected::Matrix{Symbol}
	# For integer intersection computations
	cache_int::MI
	# For floating point intersection computations
	cache_float::MF
	# color
	color::MF

end

"""
Get default values for the metadata of intersections of
a certain shape for when there is no intersection.
"""
function Intersection(
	n_intersections::Int;
	size_cache_int::Int = 1,
	size_cache_float::Int = 12,
	matrix_prototype::AbstractMatrix{F} where {F <: AbstractFloat} = zeros(Float32, 3, 3),
)::Intersection
	t = similar(matrix_prototype, (n_intersections, 1))
	t .= Inf
	return Intersection(
		Ray(n_intersections; matrix_prototype), # ray_camera
		Ray(n_intersections; matrix_prototype), # ray
		Ray(n_intersections; matrix_prototype), # ray_transformed
		t,
		fill(:none, (n_intersections, 1)), # name_intersected
		similar(matrix_prototype, Int, (n_intersections, size_cache_int)), # int_metadata_int
		similar(matrix_prototype, (n_intersections, size_cache_float)), # int_metadata_float
		similar(matrix_prototype, (n_intersections, 3)), # color
	)
end

"""
Get views of the arrays in an intersection object used for
the computations for one ray.
"""
function get_caches(intersections::Intersection, index::Int)
	ray_loc = view(intersections.ray.loc, index, :)
	ray_dir = view(intersections.ray.dir, index, :)
	ray_camera_loc = view(intersections.ray_camera.loc, index, :)
	ray_camera_dir = view(intersections.ray_camera.dir, index, :)
	t = view(intersections.t, index, :)
	cache_int = view(intersections.cache_int, index, :)
	cache_float = view(intersections.cache_float, index, :)
	color = view(intersections.color, index, :)
	name_intersected = view(intersections.name_intersected, index, :)
	return ray_loc,
	ray_dir,
	ray_camera_loc,
	ray_camera_dir,
	t,
	cache_int,
	cache_float,
	color,
	name_intersected
end

"""
Set the data to a non-intersected state.
"""
function reset_intersection!(
	t::AbstractVector,
	name_intersected::AbstractVector{Symbol},
)::Nothing
	t[1] = Inf
	name_intersected[1] = :none
	return nothing
end

"""
Compute the intersections of a ray with a sphere.
Returns (nothing, nothing) if the intersections do not exist.
"""
function intersect_sphere(
	ray_loc::VF,
	ray_dir::VF,
	Rsq::F,
)::Tuple{Union{F, Nothing}, Union{F, Nothing}} where {F <: AbstractFloat, VF <: AbstractVector{F}}
	a = norm(ray_dir)^2
	b = 2 * dot(ray_dir, ray_loc)
	c = norm(ray_loc)^2 - Rsq
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
function _intersect_ray!(
	t::AbstractVector{F},
	cache_int::AbstractVector{Int},
	cache_float::AbstractVector{F},
	ray_loc::AbstractVector{F},
	ray_dir::AbstractVector{F},
	sphere::Sphere{F},
)::Bool where {F}
	(; Rsq) = sphere
	t_int_candidate = intersect_sphere(ray_loc, ray_dir, Rsq)[1]

	closer_intersection_found = false

	if !isnothing(t_int_candidate)
		if t_int_candidate < t[1]
			closer_intersection_found = true
			t[1] = t_int_candidate
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
	t::AbstractVector{F},
	cache_int::AbstractVector{Int},
	cache_float::AbstractVector{F},
	ray_loc::AbstractVector{F},
	ray_dir::AbstractVector{F},
	fractal_shape::FractalShape{F, S, T};
	current_depth::Int = 1,
	vec_temp::Union{AbstractVector{F}, Nothing} = nothing,
)::Bool where {F, S, T}
	(; depth, shape, subshape_transforms) = fractal_shape
	closer_intersection_found = false
	t[1] = Inf

	if isnothing(vec_temp)
		vec_temp = view(cache_float, 1:3)
	end

	for subshape_transform ∈ subshape_transforms
		t[1] /= subshape_transform.scaling
		inverse_transform!(ray_loc, ray_dir, vec_temp, subshape_transform)
		if _intersect_ray!(t, cache_int, cache_float, ray_loc, ray_dir, shape)
			if current_depth < depth
				if _intersect_ray!(
					t,
					cache_int,
					cache_float,
					ray_loc,
					ray_dir,
					fractal_shape;
					current_depth = current_depth + 1,
					vec_temp,
				)
					closer_intersection_found = true
				end
			else
				closer_intersection_found = true
			end
		end
		forward_transform!(ray_loc, ray_dir, vec_temp, subshape_transform)
		t[1] *= subshape_transform.scaling
	end
	return closer_intersection_found
end

"""
Compute the intersection of a ray with a triangle given by
the triangle vertices.
"""
function _intersect_ray!(
	ray_loc::VF,
	ray_dir::VF,
	cache_float::VF,
	triangle_vertices::AbstractArray{F};
	normal::Union{AbstractVector{F}, Nothing} = nothing,
	t_int_prev::Union{F, Nothing} = nothing,
)::F where {F <: AbstractFloat, VF <: AbstractVector{F}}

	if isnothing(normal)
		u = view(cache_float, 1:3)
		v = view(cache_float, 4:6)

		c = view(triangle_vertices, 3, :)
		u = view(triangle_vertices, 1, :) - c
		v = view(triangle_vertices, 2, :) - c
		normal = cross(u, v)
		normalize!(normal)
	end

	diff = view(cache_float, 7:9)
	diff .= 0.0
	@views(diff[:] .+= triangle_vertices[3, :])
	diff .-= ray_loc
	t_int_candidate = dot(diff, normal) / dot(ray_dir, normal)

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
		u = view(cache_float, 1:3)
		v = view(cache_float, 4:6)

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
	diff2 = view(cache_float, 10:12)
	diff2 .= @. t_int_candidate * ray_dir - diff
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
	t::AbstractVector{F},
	cache_int::AbstractVector{Int},
	cache_float::VF,
	ray_loc::VF,
	ray_dir::VF,
	shape::TriangleShape{F},
)::Bool where {F <: AbstractFloat, VF <: AbstractVector{F}}
	(; vertices, faces, normals) = shape

	closer_intersection_found = false

	for i ∈ 1:shape.n_faces
		triangle_vertices = view(vertices, view(faces, i, :), :)
		normal = view(normals, i, :)
		t_int_candidate = _intersect_ray!(
			ray_loc,
			ray_dir,
			cache_float,
			triangle_vertices;
			normal,
			t_int_prev = t[1],
		)
		if t_int_candidate < t[1]
			closer_intersection_found = true
			t[1] = t_int_candidate
			cache_int[1] = i
		end
	end

	return closer_intersection_found
end

"""
Compute the gradient of a scalar field with finite differences.
"""
function ∇f_finitediff!(
	grad,
	loc,
	loc_perturbed,
	f::ScalarField;
	eps::F = 1e-4,
)::Nothing where {F <: AbstractFloat}
	loc_perturbed .= loc
	for i ∈ 1:3
		loc_perturbed[i] += eps
		grad[i] = (f(loc_perturbed) - f(loc)) / eps
		loc_perturbed[i] -= eps
	end
	return nothing
end

"""
Compute the intersection between a ray and an implicit surface.
"""
function _intersect_ray!(
	t::AbstractVector{F},
	cache_int::AbstractVector{Int},
	cache_float::VF,
	ray_loc::VF,
	ray_dir::VF,
	shape::ImplicitSurface{F},
)::Bool where {F <: AbstractFloat, VF <: AbstractVector{F}}
	(; f, ∇f!, itermax, tol, R_bound, n_divisions) = shape

	# Compute intersections of the ray with the bounding sphere
	bound_lower, bound_upper = intersect_sphere(ray_loc, ray_dir, R_bound^2)

	# If there are no intersections with the bounding sphere,
	# there are certainly no intersections with the implicit surface
	if isnothing(bound_lower)
		return false
	end

	# Compute the value of f at the closer bounding sphere intersection
	loc_int = view(cache_float, 1:3)
	loc_int .= view(ray_dir, :)
	loc_int .*= bound_lower
	loc_int .+= view(ray_loc, :)
	fval_lower = f(loc_int)

	# Find a sign change in f between the closer and further
	# bounding sphere intersections starting from the closer
	# and stepping with Δt
	t_0 = zero(F)
	Δt = (bound_upper - bound_lower) / n_divisions
	found_t_0 = false
	for i ∈ 1:n_divisions
		t_0 = bound_lower + i * Δt
		loc_int .= view(ray_dir, :)
		loc_int .*= t_0
		loc_int .+= view(ray_loc, :)
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

	grad = view(cache_float, 4:6)
	loc_perturbed = view(cache_float, 7:9)

	# Find the zero of f along the ray with the desired
	# tolerance using Newton iterations
	t_n = t_0
	for _ ∈ 1:itermax
		loc_int .= view(ray_dir, :)
		loc_int .*= t_n
		loc_int .+= view(ray_loc, :)
		fval = f(loc_int)
		error = abs(fval)
		if error < tol
			if t_n < t[1]
				t[1] = t_n
				return true
			else
				return false
			end
		else
			# If no analytical gradient of f is provided,
			# compute a finite difference gradient.
			if isnothing(∇f!)
				∇f_finitediff!(grad, loc_int, loc_perturbed, f)
			else
				∇f!(grad, loc_int)
			end
			t_n -= fval / dot(grad, ray_dir)
		end
	end
	return false
end

function _intersect_ray!(
	t::AbstractVector{F},
	cache_int::AbstractVector{Int},
	cache_float::VF,
	ray_loc::VF,
	ray_dir::VF,
	shape::RevolutionSurface{F},
)::Bool where {F <: AbstractFloat, VF <: AbstractVector{F}}
	(; z_min, z_max, r, r_max, itermax, tol, n_divisions, dr) = shape

	closer_intersection_found = false

	if ray_dir[3] ≈ zero(F)
		return closer_intersection_found
	end

	# Top disk
	t_top = (z_max - ray_loc[3]) / ray_dir[3]
	r_top = sqrt((ray_loc[1] + ray_dir[1] * t_top)^2 + (ray_loc[2] + ray_dir[2] * t_top)^2)
	if r_top <= r(z_max)
		if t_top < t[1]
			closer_intersection_found = true
			t[1] = t_top
		end
	end

	# Bottom disk
	t_bottom = (z_min - ray_loc[3]) / ray_dir[3]
	r_bottom = sqrt(
		(ray_loc[1] + ray_dir[1] * t_bottom)^2 + (ray_loc[2] + ray_dir[2] * t_bottom)^2,
	)
	if r_bottom <= r(z_min)
		if t_bottom < t[1]
			closer_intersection_found = true
			t[1] = t_bottom
		end
	end

	# Bounding cylinder
	a = ray_dir[1]^2 + ray_dir[2]^2
	b = 2 * (ray_dir[1] * ray_loc[1] + ray_dir[2] * ray_loc[2])
	c = ray_loc[1]^2 + ray_loc[2]^2
	discr = b^2 - 4 * a * (c - r_max^2)
	if discr < zero(F)
		return closer_intersection_found
	end

	denom = 2 * a
	sqrt_discr = sqrt(discr)
	t_int_min = (-b - sqrt_discr) / denom
	t_int_max = (-b + sqrt_discr) / denom

	if ray_dir[3] > 0
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
		fval = r(ray_loc[3] + ray_dir[3] * t_0) - sqrt(a * t_0^2 + b * t_0 + c)
		if fval >= 0
			found_t_0 = true
			break
		end
	end

	if !found_t_0
		return closer_intersection_found
	end

	loc_int = view(cache_float, 1:3)

	t_n = t_0
	for i ∈ 1:itermax
		loc_int .= view(ray_dir, :)
		loc_int *= t_n
		loc_int += view(ray_loc, :)
		dist = sqrt(a * t_n^2 + b * t_n + c)
		z_n = ray_loc[3] + ray_dir[3] * t_n
		fval = r(z_n) - dist
		error = abs(fval)
		if error < tol
			if t_n < t[1]
				t[1] = t_n
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
			df = ray_dir[3] * dr_value - (b / 2 + a * t_n) / dist
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
