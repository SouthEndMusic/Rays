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
