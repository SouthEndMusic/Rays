using Test
using Rays: Rays
using Random: seed!
using LinearAlgebra: normalize, norm

@testset "Bounding_box_transform" begin
    bounding_box = Rays.BoundingBox(fill(-1.0, 3), fill(1.0, 3))
    transform = Rays.translation(ones(3)) ∘ Rays.rotation(Float64[1, 0, 0], π / 4)
    bounding_box = Rays.transform_bounding_box(bounding_box, transform)
    @test bounding_box.coordinates_min ≈ [0, 1 - sqrt(2), 1 - sqrt(2)]
    @test bounding_box.coordinates_max ≈ [2, 1 + sqrt(2), 1 + sqrt(2)]
end

@testset "Bounding_box_intersection" begin
    nothing
end

@testset "Scene_partition" begin
    seed!(314156)
    n_cubes = 250
    scene = Rays.Scene()

    for i ∈ 1:n_cubes
        center = Float32.(1.2 * (rand(3) * 2 .- 1))
        R = 0.25f0 / ((2 * norm(center))^2 + 1)
        transform =
            Rays.translation(center) ∘
            Rays.rotation(normalize(rand(Float32, 3)), Float32(2π) * rand(Float32))
        cube = Rays.Cube(R)
        push!(scene, cube; transform)
    end

    Rays.partition!(scene)
    function overlap(box_1::Rays.BoundingBox, box_2::Rays.BoundingBox)::Bool
        for dim ∈ 1:3
            min_of_max = min(box_1.coordinates_max[dim], box_2.coordinates_max[dim])
            max_of_min = max(box_1.coordinates_min[dim], box_2.coordinates_min[dim])
            if max_of_min > min_of_max
                return false
            end
        end
        return true
    end

    for node ∈ scene.partition
        @test all(
            overlap(node.bounding_box, Rays.get_bounding_box(scene, name)) for
            name ∈ node.identifiers
        )
    end

    for node ∈ scene.partition
        bounding_box = node.bounding_box
        @test all(bounding_box.coordinates_max .>= bounding_box.coordinates_min)
    end

    names_in_endnodes = Set{Symbol}()
    for node ∈ scene.partition
        if isempty(node.child_indices)
            union!(names_in_endnodes, node.identifiers)
        end
    end
    @test names_in_endnodes == Set(keys(scene.shapes))

    for node ∈ scene.partition
        if node.dim_split != 0
            child_node_1 = scene.partition[node.child_indices[1]]
            child_node_2 = scene.partition[node.child_indices[2]]
            @test child_node_1.bounding_box.coordinates_max[node.dim_split] ==
                  child_node_2.bounding_box.coordinates_min[node.dim_split]
        end
    end
end

@testset "TriangleShape_partition" begin
    seed!(1415)

    scene = Rays.Scene()
    camera = Rays.Camera(; screen_res = (100, 100))
    push!(scene, camera)
    from = Float32[2.5, 0, 0]
    to = zeros(Float32, 3)
    Rays.look_at!(camera, from, to)

    N = 100
    vertices = zeros(Float32, 3 * N, 3)
    faces = zeros(Int, N, 3)
    face_loc = zeros(Float32, 3)
    for i ∈ 1:N
        offset = (i - 1) * 3
        faces[i, :] = [1, 2, 3] .+ offset
        face_loc[2:3] = 2.0f0 * rand(Float32, 2) .- 1.0f0
        for j ∈ 1:3
            perturbation = 0.1 * (2.0f0 * rand(Float32, 3) .- 1.0f0)
            vertices[offset+j, :] .= face_loc + perturbation
        end
    end
    triangle_shape = Rays.TriangleShape(vertices, faces, name = :donut)
    push!(scene, triangle_shape)

    @test isempty(triangle_shape.partition)
    Rays.render!(scene)
    canvas_without_partition = copy(camera.canvas[1, :, :])

    Rays.partition!(triangle_shape)
    @test !isempty(triangle_shape.partition)
    Rays.render!(scene)
    canvas_with_partition = camera.canvas[1, :, :]
    @test canvas_without_partition == canvas_with_partition
end
