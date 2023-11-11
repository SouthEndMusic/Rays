using Test
using Rays: Rays
using Random: seed!

@testset "blur_kernel" begin
    kernel, Δi_max = Rays.get_blur_kernel(1.6)
    @test sum(kernel) ≈ 1.0
    @test all(kernel .>= 0.0)
    @test kernel ≈ Float32[0.00029104948, 0.2252549, 0.5489081, 0.2252549, 0.00029104948]
end

@testset "shape_renders" begin
    scene = Rays.Scene()

    camera = Rays.Camera(; screen_res = [25, 25])
    from = Float32[2.0, 2.0, 2.0]
    to = zeros(Float32, 3)
    Rays.look_at!(camera, from, to)
    push!(scene, camera)

    origin = zeros(Float32, 3)
    R = 1.0f0

    function my_field(loc::Vector{F})::F where {F}
        x, y, z = loc
        out = -4.0f0
        for i ∈ 0:3
            θ = convert(F, π * (i / 2 + 0.33))
            out += one(F) / sqrt((x - cos(θ))^2 + (y - sin(θ))^2 + z^2)
        end
        return out
    end

    shapes = [
        Rays.Cube(origin, R),
        Rays.menger_sponge(origin, R, 4),
        Rays.Sphere(origin, R),
        Rays.Tetrahedron(origin, R),
        Rays.sierpinski_pyramid(origin, R, 4),
        Rays.ImplicitSurface(
            my_field,
            origin,
            R_bound = 1.5f0,
            n_divisions = 50,
            tol = 1.0f-5,
            itermax = 10,
            name = :equipotential_surface,
        ),
    ]

    push!(scene, shapes[1])
    Rays.set_dropoff_curve_default!(scene, camera)

    function simple_view!(shape)
        Rays.clear_shapes!(scene)
        push!(scene, shape)
        Rays.render!(scene)
        return nothing
    end

    simple_view!(shapes[1])
    @test camera.canvas[1, :, :] ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.26033032 0.2620107 0.2625717 0.2620107 0.26033032 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.33769637 0.34310693 0.3475672 0.35105753 0.35356224 0.35506976 0.355573 0.35506976 0.35356224 0.35105753 0.3475672 0.34310693 0.33769637 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.3850698 0.39318103 0.40053165 0.40709293 0.4128387 0.41774482 0.42178982 0.4249556 0.42722774 0.42859524 0.42905188 0.42859524 0.42722774 0.4249556 0.42178982 0.41774482 0.4128387 0.40709293 0.40053165 0.39318103 0.3850698 0.0 0.0
        0.0 0.4107151 0.44810504 0.45551854 0.46223855 0.4682383 0.47349334 0.47798127 0.48168206 0.48457885 0.48665816 0.48790967 0.4883275 0.48790967 0.48665816 0.48457885 0.48168206 0.47798127 0.47349334 0.4682383 0.46223855 0.45551854 0.44810504 0.4107151 0.0
        0.0 0.3846643 0.44395256 0.4937718 0.5129163 0.518441 0.52328086 0.5274149 0.5308243 0.5334934 0.53540933 0.53656256 0.5369476 0.53656256 0.53540933 0.5334934 0.5308243 0.5274149 0.52328086 0.518441 0.5129163 0.4937718 0.44395256 0.3846643 0.0
        0.0 0.0 0.41971773 0.47369426 0.51923066 0.5580049 0.56470406 0.5685332 0.5716915 0.5741643 0.5759394 0.577008 0.5773648 0.577008 0.5759394 0.5741643 0.5716915 0.5685332 0.56470406 0.5580049 0.51923066 0.47369426 0.41971773 0.0 0.0
        0.0 0.0 0.3921659 0.45095003 0.50020325 0.54189885 0.5775019 0.60311294 0.6060517 0.60835254 0.6100044 0.610999 0.611331 0.610999 0.6100044 0.60835254 0.6060517 0.60311294 0.5775019 0.54189885 0.50020325 0.45095003 0.3921659 0.0 0.0
        0.0 0.0 0.3607419 0.42512685 0.47866756 0.5237062 0.5619571 0.5946996 0.62290955 0.6373553 0.6388979 0.63982654 0.64013666 0.63982654 0.6388979 0.6373553 0.62290955 0.5946996 0.5619571 0.5237062 0.47866756 0.42512685 0.3607419 0.0 0.0
        0.0 0.0 0.32476962 0.3957315 0.4542519 0.5031387 0.54441583 0.5795733 0.6097352 0.63576627 0.6583444 0.66445875 0.66474915 0.66445875 0.6583444 0.63576627 0.6097352 0.5795733 0.54441583 0.5031387 0.4542519 0.3957315 0.32476962 0.0 0.0
        0.0 0.0 0.0 0.36217082 0.4265167 0.47986096 0.5246147 0.56252867 0.5949069 0.62274206 0.6468043 0.6676986 0.68591166 0.6676986 0.6468043 0.62274206 0.5949069 0.56252867 0.5246147 0.47986096 0.4265167 0.36217082 0.0 0.0 0.0
        0.0 0.0 0.0 0.3237158 0.3949327 0.45347494 0.5022464 0.5433225 0.5782285 0.60811114 0.6338519 0.6561363 0.67551124 0.6561363 0.6338519 0.60811114 0.5782285 0.5433225 0.5022464 0.45347494 0.3949327 0.3237158 0.0 0.0 0.0
        0.0 0.0 0.0 0.27946335 0.35885668 0.42350644 0.4769504 0.5216736 0.5594752 0.5916912 0.61933637 0.64319277 0.6638788 0.64319277 0.61933637 0.5916912 0.5594752 0.5216736 0.4769504 0.42350644 0.35885668 0.27946335 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.31750453 0.38938814 0.44830257 0.49725574 0.5383905 0.5732758 0.603089 0.62872815 0.65089625 0.62872815 0.603089 0.5732758 0.5383905 0.49725574 0.44830257 0.38938814 0.31750453 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.26989847 0.35042745 0.4157945 0.46968466 0.5146762 0.5526285 0.58491886 0.61258554 0.6364333 0.61258554 0.58491886 0.5526285 0.5146762 0.46968466 0.4157945 0.35042745 0.26989847 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.21480328 0.30576903 0.3788103 0.43850183 0.48798168 0.5294751 0.5646068 0.59458697 0.62034297 0.59458697 0.5646068 0.5294751 0.48798168 0.43850183 0.3788103 0.30576903 0.21480328 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.25435448 0.3366047 0.4031638 0.45789844 0.5035005 0.5419047 0.57453316 0.6024629 0.57453316 0.5419047 0.5035005 0.45789844 0.4031638 0.3366047 0.25435448 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.28825128 0.36300737 0.42393595 0.47433197 0.51652277 0.5521945 0.58260745 0.5521945 0.51652277 0.47433197 0.42393595 0.36300737 0.28825128 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.31721908 0.38550377 0.44152784 0.48812193 0.5273049 0.5605646 0.5273049 0.48812193 0.44152784 0.38550377 0.31721908 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.34188253 0.40455896 0.45630252 0.49955547 0.5360903 0.49955547 0.45630252 0.40455896 0.34188253 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.36279106 0.4205929 0.4685865 0.5089046 0.4685865 0.4205929 0.36279106 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3804204 0.43396723 0.4786756 0.43396723 0.3804204 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.39518422 0.44501364 0.39518422 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    simple_view!(shapes[2])
    @test camera.canvas[1, :, :] ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.25878143 0.2620108 0.0 0.2620108 0.25878143 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.33559555 0.34148628 0.33424973 0.3304587 0.35354263 0.34970015 0.35557306 0.34970015 0.35354263 0.3304587 0.33424973 0.34148628 0.33559555 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.38506985 0.39318103 0.3937992 0.407093 0.41283882 0.41774487 0.42178988 0.41989815 0.39886725 0.37514722 0.3483848 0.37514722 0.39886725 0.41989815 0.42178988 0.41774487 0.41283882 0.407093 0.3937992 0.39318103 0.38506985 0.0 0.0
        0.0 0.4107151 0.44810504 0.45551556 0.46112627 0.46657795 0.4734934 0.47798133 0.48168212 0.4845789 0.48665816 0.3663556 0.0 0.3663556 0.48665816 0.4845789 0.48168212 0.47798133 0.4734934 0.46657795 0.46112627 0.45551556 0.44810504 0.4107151 0.0
        0.0 0.38466442 0.4433506 0.49377185 0.5129163 0.5184411 0.5232809 0.4969778 0.49157983 0.52995515 0.5314938 0.53230137 0.5323664 0.53230137 0.5314938 0.52995515 0.49157983 0.4969778 0.5232809 0.5184411 0.5129163 0.49377185 0.4433506 0.38466442 0.0
        0.0 0.0 0.39776057 0.47369432 0.5192307 0.55800486 0.5647041 0.5639129 0.57169163 0.5741644 0.57593954 0.577008 0.5773649 0.577008 0.57593954 0.5741644 0.57169163 0.5639129 0.5647041 0.55800486 0.5192307 0.47369432 0.39776057 0.0 0.0
        0.0 0.0 0.38938498 0.4459986 0.5002033 0.5324459 0.5775019 0.60311294 0.6001172 0.6083526 0.61000454 0.5931915 0.5735719 0.5931915 0.61000454 0.6083526 0.6001172 0.60311294 0.5775019 0.5324459 0.5002033 0.4459986 0.38938498 0.0 0.0
        0.0 0.0 0.3607419 0.42512685 0.47760534 0.5008365 0.5619572 0.5946997 0.6185523 0.6373553 0.6232054 0.6398266 0.64013666 0.6398266 0.6232054 0.6373553 0.6185523 0.5946997 0.5619572 0.5008365 0.47760534 0.42512685 0.3607419 0.0 0.0
        0.0 0.0 0.32476962 0.39573157 0.45425195 0.4910217 0.5444159 0.5795733 0.60973525 0.63576627 0.6583444 0.6644588 0.64717627 0.6644588 0.6583444 0.63576627 0.60973525 0.5795733 0.5444159 0.4910217 0.45425195 0.39573157 0.32476962 0.0 0.0
        0.0 0.0 0.0 0.36217093 0.4265167 0.47986108 0.5246147 0.5625287 0.5949069 0.6227422 0.64680433 0.6676986 0.68591166 0.6676986 0.64680433 0.6227422 0.5949069 0.5625287 0.5246147 0.47986108 0.4265167 0.36217093 0.0 0.0 0.0
        0.0 0.0 0.0 0.32371587 0.3949327 0.4093827 0.5022465 0.5425695 0.5782286 0.5872028 0.6338519 0.65613633 0.67551124 0.65613633 0.6338519 0.5872028 0.5782286 0.5425695 0.5022465 0.4093827 0.3949327 0.32371587 0.0 0.0 0.0
        0.0 0.0 0.0 0.2794633 0.35885674 0.3940444 0.37502414 0.5216737 0.5594753 0.583907 0.5784467 0.6402279 0.6638788 0.6402279 0.5784467 0.583907 0.5594753 0.5216737 0.37502414 0.3940444 0.35885674 0.2794633 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.31750453 0.37707967 0.35682774 0.33963847 0.53839064 0.5685344 0.603089 0.6278858 0.6508963 0.6278858 0.603089 0.5685344 0.53839064 0.33963847 0.35682774 0.37707967 0.31750453 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.25477654 0.3504274 0.35837495 0.36454868 0.5146763 0.5526285 0.58491886 0.6123476 0.6364333 0.6123476 0.58491886 0.5526285 0.5146763 0.36454868 0.35837495 0.3504274 0.25477654 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.21480334 0.30576897 0.37881035 0.3810994 0.48389983 0.5294752 0.5646068 0.5945871 0.62034297 0.5945871 0.5646068 0.5294752 0.48389983 0.3810994 0.37881035 0.30576897 0.21480334 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.2543546 0.3366047 0.40316385 0.45789844 0.5035006 0.54190475 0.56824887 0.6024629 0.56824887 0.54190475 0.5035006 0.45789844 0.40316385 0.3366047 0.2543546 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.28717977 0.32407427 0.42393607 0.47433197 0.496857 0.5521946 0.58260745 0.5521946 0.496857 0.47433197 0.42393607 0.32407427 0.28717977 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.30846572 0.38550377 0.4415279 0.47916585 0.527305 0.56056464 0.527305 0.47916585 0.4415279 0.38550377 0.30846572 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.33190733 0.38743925 0.44159567 0.49955547 0.5360904 0.49955547 0.44159567 0.38743925 0.33190733 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.35174626 0.3878783 0.46858656 0.50890464 0.46858656 0.3878783 0.35174626 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.37349015 0.43396723 0.47867566 0.43396723 0.37349015 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.39518428 0.44501364 0.39518428 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    simple_view!(shapes[3])
    @test camera.canvas[1, :, :] ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.4197964 0.42687243 0.4197964 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3969791 0.44742823 0.4656698 0.47466773 0.47744614 0.47466773 0.4656698 0.44742823 0.3969791 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.41979516 0.4624039 0.4827295 0.49466258 0.50113606 0.5032002 0.50113606 0.49466258 0.4827295 0.4624039 0.41979516 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.39698702 0.46240497 0.4876992 0.50320035 0.5129284 0.51835704 0.52010775 0.51835704 0.5129284 0.50320035 0.4876992 0.46240497 0.39698702 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.44742912 0.4827295 0.5032002 0.51657724 0.5252001 0.5300751 0.53165585 0.5300751 0.5252001 0.51657724 0.5032002 0.4827295 0.44742912 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.46567047 0.49466228 0.51292765 0.5252001 0.533216 0.5377786 0.5392625 0.5377786 0.533216 0.5252001 0.51292765 0.49466228 0.46567047 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.41979462 0.4746678 0.50113606 0.5183563 0.5300748 0.53777874 0.54217833 0.5436114 0.54217833 0.53777874 0.5300748 0.5183563 0.50113606 0.4746678 0.41979462 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.42687172 0.4774462 0.50320035 0.52010727 0.53165567 0.53926253 0.5436114 0.54502845 0.5436114 0.53926253 0.53165567 0.52010727 0.50320035 0.4774462 0.42687172 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.4197964 0.4746685 0.50113606 0.51835644 0.53007513 0.5377786 0.5421783 0.54361117 0.5421783 0.5377786 0.53007513 0.51835644 0.50113606 0.4746685 0.4197964 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.46566993 0.49466228 0.51292765 0.5252001 0.533216 0.5377785 0.53926253 0.5377785 0.533216 0.5252001 0.51292765 0.49466228 0.46566993 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.4474296 0.4827295 0.5032002 0.5165772 0.52520037 0.53007495 0.53165585 0.53007495 0.52520037 0.5165772 0.5032002 0.4827295 0.4474296 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.39698702 0.46240497 0.48769963 0.50320035 0.5129282 0.51835704 0.52010787 0.51835704 0.5129282 0.50320035 0.48769963 0.46240497 0.39698702 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.41979218 0.4624039 0.4827289 0.4946627 0.50113606 0.5032002 0.50113606 0.4946627 0.4827289 0.4624039 0.41979218 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.39697808 0.44742823 0.4656698 0.47466755 0.47744578 0.47466755 0.4656698 0.44742823 0.39697808 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.4197964 0.42687243 0.4197964 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    simple_view!(shapes[4])
    @test camera.canvas[1, :, :] ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.44133765 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.43885893 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.4411173 0.435412 0.42861986 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.43676084 0.43093586 0.42400026 0.4159186 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.4313305 0.4253695 0.41827494 0.41001004 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.42976797 0.4247682 0.41865373 0.41138315 0.4029187 0.3932252 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.42216963 0.41701776 0.41073084 0.40326607 0.39458424 0.38464934 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.41335 0.4080249 0.4015451 0.39386648 0.38494837 0.3747536 0.3632493 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.407654 0.40325844 0.39773762 0.391043 0.38312918 0.3739543 0.3634796 0.3516711 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.39644533 0.39184707 0.38610643 0.37917352 0.3710019 0.3615482 0.35077208 0.33863813 0.32511282 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.38389742 0.3790704 0.37308425 0.36588818 0.35743415 0.34767777 0.33657682 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3739298 0.3699674 0.36488378 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    simple_view!(shapes[5])
    @test camera.canvas[1, :, :] ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.44133765 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.4076379 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.44111723 0.0 0.41273767 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.43676084 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.43133044 0.36151457 0.41827494 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.41392207 0.0 0.41865373 0.0 0.4029187 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.27081424 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.41334993 0.27659273 0.26872313 0.25938994 0.0 0.0 0.35615575 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.407654 0.0 0.39773774 0.0 0.0 0.3739543 0.27532375 0.31934834 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3294558 0.3243721 0.0 0.0 0.0 0.0 0.0 0.0 0.32511288 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.38389748 0.0 0.33833462 0.36588818 0.33963877 0.0 0.33657682 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.35656726 0.3699674 0.36488378 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    simple_view!(shapes[6])
    @test camera.canvas[1, :, :] ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.23853612 0.2673036 0.27096063 0.2563753 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.28184754 0.2978757 0.30076307 0.29349506 0.27058768 0.0 0.0 0.0 0.0 0.30799997 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.26329798 0.3063473 0.31757885 0.31921834 0.31450444 0.3044532 0.2947163 0.32950854 0.3550362 0.36570942 0.36676592 0.35870326 0.33573985 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.4016654 0.41557038 0.41721642 0.40935516 0.38802528 0.3490113 0.33758318 0.3336929 0.32969362 0.32922477 0.34614414 0.36889172 0.382316 0.387614 0.3863358 0.37850893 0.361753 0.31803572 0.0 0.0 0.0
        0.0 0.0 0.0 0.41609913 0.43416858 0.44209838 0.4437036 0.439763 0.42976785 0.41140676 0.37906772 0.35017014 0.3420171 0.3513131 0.37535572 0.39210165 0.3993562 0.40020734 0.39613372 0.38689703 0.3702708 0.3358528 0.0 0.0 0.0
        0.0 0.0 0.40218443 0.43319786 0.44726765 0.45424694 0.45625597 0.4540751 0.44807172 0.4385503 0.42584205 0.40943164 0.35257828 0.47210592 0.47900242 0.4486314 0.4191484 0.4085381 0.3991139 0.38686508 0.36782277 0.32643032 0.0 0.0 0.0
        0.0 0.0 0.41191626 0.43840635 0.45165431 0.4585594 0.46103883 0.4601171 0.4570105 0.45466554 0.46312207 0.4954096 0.51999617 0.53177685 0.534965 0.53027 0.51409984 0.42066395 0.39379603 0.37619138 0.35032398 0.0 0.0 0.0 0.0
        0.0 0.0 0.40543234 0.43492538 0.44888556 0.45620406 0.4591546 0.45916677 0.45857733 0.46484035 0.4961183 0.5274288 0.5446235 0.5533471 0.5559839 0.5530468 0.543476 0.5219152 0.36545295 0.33643103 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.42038542 0.43766493 0.4460765 0.44938827 0.44957966 0.45028865 0.46833438 0.51571333 0.5426761 0.55726886 0.5648787 0.56731397 0.5650224 0.5574244 0.54195714 0.50116205 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.40888858 0.42183816 0.4250151 0.42113864 0.0 0.0 0.52397114 0.5495378 0.56324875 0.57047147 0.57287663 0.5709073 0.5641493 0.5507455 0.522255 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5228853 0.5501761 0.5641255 0.57144296 0.5739565 0.5721606 0.5656954 0.5528492 0.5261139 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.50765264 0.54452586 0.56008625 0.5680294 0.5708202 0.56909585 0.56248105 0.5489912 0.5179303 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.52925634 0.5499463 0.55950165 0.5628723 0.56110084 0.55368066 0.5373781 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.52783704 0.54284614 0.5476992 0.54558456 0.5352979 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.50652003 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    # Many shapes
    seed!(31415)
    Rays.clear_shapes!(scene)
    n_cubes = 250

    for i = 1:n_cubes
        center = rand(Float32, 3) * 2 .- 1
        R = rand(Float32) / 10
        cube = Rays.Cube(center, R)
        push!(scene, cube)
    end

    Rays.render!(scene; name_camera = :camera)
    @test camera.canvas[1, :, :] ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.23247516 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.22392458 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.35961974 0.35940945 0.0 0.42090762 0.0 0.0 0.0 0.0 0.0 0.38185865 0.0 0.41715014 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.36851722 0.35273343 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.30387497 0.0 0.0 0.0 0.0 0.0 0.0 0.40336895 0.0 0.0
        0.0 0.0 0.0 0.46761554 0.47913545 0.0 0.49140793 0.3655681 0.0 0.46505684 0.0 0.0 0.39685154 0.0 0.2915817 0.40757644 0.2810868 0.0 0.34895468 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.45094156 0.47547925 0.45843506 0.47332096 0.5344341 0.0 0.45625287 0.26728523 0.27278107 0.43932432 0.40785366 0.2770154 0.27487588 0.25882626 0.0 0.36548775 0.498155 0.49231714 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.47768325 0.46611488 0.42576474 0.49111634 0.56150806 0.50120693 0.506639 0.0 0.0 0.0 0.4007219 0.0 0.20124161 0.0 0.48115796 0.375665 0.5208099 0.48406833 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.4691493 0.0 0.44843078 0.4840421 0.50133634 0.51063657 0.5331575 0.5138581 0.0 0.2297579 0.0 0.6080628 0.60638857 0.0 0.47126138 0.4888779 0.51714253 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.39501137 0.39789546 0.46933913 0.43442768 0.38242626 0.5210206 0.0 0.40485168 0.6649281 0.66463786 0.6305157 0.6218624 0.2913131 0.0 0.49472708 0.4847741 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.38171715 0.3856457 0.0 0.0 0.57724035 0.5914938 0.5888876 0.0 0.6810024 0.6858075 0.67260563 0.6083668 0.2633593 0.40163833 0.35089302 0.30773014 0.4387691 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.3260072 0.3759023 0.47794068 0.0 0.0 0.0 0.5909981 0.5767 0.5275949 0.6704486 0.68128896 0.66055155 0.44995713 0.45496947 0.42659473 0.4119832 0.29342854 0.3122275 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.327999 0.37552333 0.37552673 0.3754937 0.3391018 0.18871856 0.0 0.41468263 0.41708636 0.27932805 0.6692454 0.45187855 0.28832614 0.43593824 0.41052407 0.41751575 0.3613782 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.32910657 0.0 0.32528013 0.35039312 0.30674314 0.35013366 0.33113706 0.4101668 0.28594398 0.37328982 0.22475725 0.0 0.3175357 0.45466095 0.38548297 0.23613697 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.366589 0.47300333 0.46554685 0.46863234 0.33452886 0.0 0.42574596 0.4252612 0.36655605 0.46677077 0.46967006 0.2387374 0.4003831 0.38197714 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.4654994 0.46631742 0.4634207 0.0 0.38763064 0.4143737 0.40698242 0.33552355 0.42802596 0.4488933 0.37284076 0.3579901 0.30623132 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.38360667 0.0 0.46178812 0.3595655 0.3576827 0.3354419 0.36104828 0.0 0.26618946 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.28747225 0.0 0.0 0.0 0.4705739 0.47963756 0.34289747 0.4122072 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.42031586 0.30041307 0.5432296 0.4969417 0.4024781 0.42992884 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.43399358 0.43612236 0.41289294 0.5398383 0.50607157 0.35461897 0.4076206 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.422495 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.38392496 0.3852446 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.37669802 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
end
