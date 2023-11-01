import Rays

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

    shapes = [
        Rays.Cube(origin, R),
        Rays.menger_sponge(origin, R, 4),
        Rays.Sphere(origin, R),
        Rays.Tetrahedron(origin, R),
        Rays.sierpinski_pyramid(origin, R, 4),
    ]

    function simple_view(shape)
        push!(scene, shape)
        Rays.shape_view!(scene)
        Rays.cam_is_source!(camera)
        delete!(scene.shapes, shape.name)
        return nothing
    end

    simple_view(shapes[1])
    @test camera.canvas_grayscale ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.09663821 0.100205116 0.10139587 0.100205116 0.09663821 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.26085952 0.27234426 0.28181195 0.2892207 0.2945373 0.2977372 0.29880545 0.2977372 0.2945373 0.2892207 0.28181195 0.27234426 0.26085952 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.3614169 0.3786342 0.39423713 0.40816444 0.42036074 0.43077463 0.4393609 0.44608068 0.45090368 0.45380634 0.45477563 0.45380634 0.45090368 0.44608068 0.4393609 0.43077463 0.42036074 0.40816444 0.39423713 0.3786342 0.3614169 0.0 0.0
        0.0 0.41585302 0.49521887 0.5109552 0.5252195 0.53795487 0.54910946 0.5586358 0.56649137 0.57264024 0.5770539 0.57971036 0.58059734 0.57971036 0.5770539 0.57264024 0.56649137 0.5586358 0.54910946 0.53795487 0.5252195 0.5109552 0.49521887 0.41585302 0.0
        0.0 0.36055622 0.48640472 0.59215367 0.6327908 0.6445179 0.6547911 0.6635663 0.67080325 0.67646885 0.6805357 0.6829836 0.68380105 0.6829836 0.6805357 0.67646885 0.67080325 0.6635663 0.6547911 0.6445179 0.6327908 0.59215367 0.48640472 0.36055622 0.0
        0.0 0.0 0.43496254 0.54953605 0.64619404 0.72849834 0.7427183 0.7508462 0.7575503 0.76279914 0.7665672 0.7688354 0.76959276 0.7688354 0.7665672 0.76279914 0.7575503 0.7508462 0.7427183 0.72849834 0.64619404 0.54953605 0.43496254 0.0 0.0
        0.0 0.0 0.37647945 0.50125784 0.60580546 0.6943108 0.76988363 0.824247 0.830485 0.83536893 0.83887535 0.8409863 0.8416911 0.8409863 0.83887535 0.83536893 0.830485 0.824247 0.76988363 0.6943108 0.60580546 0.50125784 0.37647945 0.0 0.0
        0.0 0.0 0.30977732 0.44644415 0.5600925 0.65569407 0.7368875 0.80638856 0.8662684 0.8969318 0.9002061 0.9021774 0.90283555 0.9021774 0.9002061 0.8969318 0.8662684 0.80638856 0.7368875 0.65569407 0.5600925 0.44644415 0.30977732 0.0 0.0
        0.0 0.0 0.23342048 0.38404804 0.5082665 0.6120364 0.6996534 0.7742805 0.83830374 0.89355874 0.94148433 0.95446306 0.95507944 0.95446306 0.94148433 0.89355874 0.83830374 0.7742805 0.6996534 0.6120364 0.5082665 0.38404804 0.23342048 0.0 0.0
        0.0 0.0 0.0 0.3128104 0.44939443 0.56262577 0.6576224 0.7381007 0.80682856 0.8659129 0.9169886 0.96134007 1.0 0.96134007 0.9169886 0.8659129 0.80682856 0.7381007 0.6576224 0.56262577 0.44939443 0.3128104 0.0 0.0 0.0
        0.0 0.0 0.0 0.23118359 0.38235238 0.5066173 0.61014235 0.6973326 0.7714259 0.83485657 0.8894952 0.9367972 0.9779236 0.9367972 0.8894952 0.83485657 0.7714259 0.6973326 0.61014235 0.5066173 0.38235238 0.23118359 0.0 0.0 0.0
        0.0 0.0 0.0 0.13725102 0.30577552 0.44300464 0.5564477 0.6513795 0.73161924 0.80000263 0.8586837 0.9093227 0.95323193 0.9093227 0.8586837 0.80000263 0.73161924 0.6513795 0.5564477 0.44300464 0.30577552 0.13725102 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.21799931 0.37058324 0.49563825 0.5995489 0.6868637 0.7609131 0.82419616 0.87861925 0.9256745 0.87861925 0.82419616 0.7609131 0.6868637 0.5995489 0.49563825 0.37058324 0.21799931 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.11694812 0.28788316 0.42663485 0.541025 0.6365264 0.717086 0.7856273 0.84435403 0.89497465 0.84435403 0.7856273 0.717086 0.6365264 0.541025 0.42663485 0.28788316 0.11694812 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.19308884 0.34813017 0.47483465 0.57986325 0.6679395 0.7425118 0.8061494 0.86082035 0.8061494 0.7425118 0.6679395 0.57986325 0.47483465 0.34813017 0.19308884 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.08395356 0.2585423 0.3998242 0.5160069 0.61280435 0.6943231 0.7635821 0.8228671 0.7635821 0.6943231 0.61280435 0.5160069 0.3998242 0.2585423 0.08395356 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.15590473 0.31458598 0.44391635 0.5508897 0.64044607 0.71616477 0.780721 0.71616477 0.64044607 0.5508897 0.44391635 0.31458598 0.15590473 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.21739331 0.36233807 0.48125786 0.5801609 0.6633329 0.7339316 0.6633329 0.5801609 0.48125786 0.36233807 0.21739331 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.26974526 0.4027857 0.5126193 0.60443044 0.6819812 0.60443044 0.5126193 0.4027857 0.26974526 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.31412688 0.43682015 0.538694 0.62427527 0.538694 0.43682015 0.31412688 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.35154787 0.46520928 0.5601097 0.46520928 0.35154787 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3828864 0.48865697 0.3828864 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    simple_view(shapes[2])
    @test camera.canvas_grayscale ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.09335021 0.10020508 0.0 0.10020508 0.09335021 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.25640002 0.2689041 0.25354332 0.24549633 0.29449555 0.2863393 0.2988054 0.2863393 0.29449555 0.24549633 0.25354332 0.2689041 0.25640002 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.36141685 0.37863418 0.3799463 0.4081644 0.4203608 0.4307747 0.43936086 0.43534535 0.39070398 0.34035456 0.28354728 0.34035456 0.39070398 0.43534535 0.43936086 0.4307747 0.4203608 0.4081644 0.3799463 0.37863418 0.36141685 0.0 0.0
        0.0 0.41585296 0.49521884 0.51094884 0.5228583 0.53443044 0.5491095 0.5586359 0.56649137 0.57264024 0.5770538 0.32169297 0.0 0.32169297 0.5770538 0.57264024 0.56649137 0.5586359 0.5491095 0.53443044 0.5228583 0.51094884 0.49521884 0.41585296 0.0
        0.0 0.36055627 0.48512685 0.59215367 0.6327908 0.64451796 0.65479124 0.59895885 0.58750075 0.6689583 0.6722243 0.6739385 0.6740766 0.6739385 0.6722243 0.6689583 0.58750075 0.59895885 0.65479124 0.64451796 0.6327908 0.59215367 0.48512685 0.36055627 0.0
        0.0 0.0 0.3883549 0.54953605 0.64619404 0.7284981 0.74271834 0.7410389 0.75755036 0.76279914 0.7665672 0.7688353 0.7695929 0.7688353 0.7665672 0.76279914 0.75755036 0.7410389 0.74271834 0.7284981 0.64619404 0.54953605 0.3883549 0.0 0.0
        0.0 0.0 0.3705765 0.49074754 0.6058054 0.6742453 0.7698836 0.82424694 0.817888 0.835369 0.8388754 0.80318725 0.7615416 0.80318725 0.8388754 0.835369 0.817888 0.82424694 0.7698836 0.6742453 0.6058054 0.49074754 0.3705765 0.0 0.0
        0.0 0.0 0.30977717 0.44644412 0.5578377 0.6071495 0.7368876 0.80638856 0.8570195 0.89693177 0.86689645 0.90217733 0.90283555 0.90217733 0.86689645 0.89693177 0.8570195 0.80638856 0.7368876 0.6071495 0.5578377 0.44644412 0.30977717 0.0 0.0
        0.0 0.0 0.23342031 0.3840481 0.50826657 0.586316 0.69965345 0.77428037 0.8383038 0.8935587 0.9414843 0.954463 0.91777813 0.954463 0.9414843 0.8935587 0.8383038 0.77428037 0.69965345 0.586316 0.50826657 0.3840481 0.23342031 0.0 0.0
        0.0 0.0 0.0 0.31281042 0.44939432 0.56262594 0.6576224 0.7381007 0.8068285 0.86591303 0.91698873 0.96133995 1.0 0.96133995 0.91698873 0.86591303 0.8068285 0.7381007 0.6576224 0.56262594 0.44939432 0.31281042 0.0 0.0 0.0
        0.0 0.0 0.0 0.23118363 0.38235235 0.41302472 0.61014235 0.69573426 0.771426 0.7904752 0.8894952 0.9367973 0.9779236 0.9367973 0.8894952 0.7904752 0.771426 0.69573426 0.61014235 0.41302472 0.38235235 0.23118363 0.0 0.0 0.0
        0.0 0.0 0.0 0.1372507 0.30577558 0.38046676 0.34009334 0.6513795 0.73161936 0.78347933 0.77188903 0.90302914 0.9532318 0.90302914 0.77188903 0.78347933 0.73161936 0.6513795 0.34009334 0.38046676 0.30577558 0.1372507 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.21799915 0.34445646 0.30146867 0.2649818 0.68686384 0.7508486 0.8241962 0.87683123 0.92567456 0.87683123 0.8241962 0.7508486 0.68686384 0.2649818 0.30146867 0.34445646 0.21799915 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.084849365 0.28788292 0.30475292 0.31785753 0.63652647 0.717086 0.78562725 0.8438489 0.8949746 0.8438489 0.78562725 0.717086 0.63652647 0.31785753 0.30475292 0.28788292 0.084849365 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 4.869759f-8 0.19308858 0.3481302 0.35298905 0.57119876 0.66793954 0.74251175 0.80614954 0.8608204 0.80614954 0.74251175 0.66793954 0.57119876 0.35298905 0.3481302 0.19308858 4.869759f-8 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.08395352 0.25854224 0.39982426 0.5160069 0.61280435 0.6943231 0.7502426 0.82286704 0.7502426 0.6943231 0.61280435 0.5160069 0.39982426 0.25854224 0.08395352 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.15363023 0.23194438 0.4439164 0.55088955 0.5987023 0.71616495 0.78072095 0.71616495 0.5987023 0.55088955 0.4439164 0.23194438 0.15363023 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.19881287 0.36233804 0.48125783 0.5611502 0.6633329 0.7339317 0.6633329 0.5611502 0.48125783 0.36233804 0.19881287 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.24857129 0.36644635 0.48140168 0.6044303 0.6819812 0.6044303 0.48140168 0.36644635 0.24857129 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.29068244 0.36737823 0.538694 0.6242753 0.538694 0.36737823 0.29068244 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.33683723 0.46520916 0.5601097 0.46520916 0.33683723 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.38288635 0.48865685 0.38288635 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    simple_view(shapes[3])
    @test camera.canvas_grayscale ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.15412536 0.2019202 0.15412536 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 7.1281447f-6 0.34076345 0.46397528 0.5247514 0.54351825 0.5247514 0.46397528 0.34076345 7.1281447f-6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.154117 0.44191584 0.57920456 0.6598059 0.70353085 0.71747255 0.70353085 0.6598059 0.57920456 0.44191584 0.154117 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 6.043427f-5 0.4419233 0.6127719 0.7174738 0.78318167 0.81984884 0.83167446 0.81984884 0.78318167 0.7174738 0.6127719 0.4419233 6.043427f-5 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.34076965 0.57920456 0.71747255 0.8078277 0.8660705 0.8989982 0.90967524 0.8989982 0.8660705 0.8078277 0.71747255 0.57920456 0.34076965 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.46397993 0.65980405 0.7831767 0.8660702 0.9202131 0.9510315 0.96105397 0.9510315 0.9202131 0.8660702 0.7831767 0.65980405 0.46397993 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.15411328 0.5247517 0.70353055 0.8198442 0.89899606 0.9510321 0.9807491 0.99042875 0.9807491 0.9510321 0.89899606 0.8198442 0.70353055 0.5247517 0.15411328 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.20191555 0.54351854 0.7174738 0.831671 0.909674 0.9610543 0.99042875 1.0 0.99042875 0.9610543 0.909674 0.831671 0.7174738 0.54351854 0.20191555 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.15412536 0.52475667 0.70353085 0.81984484 0.8989985 0.9510315 0.9807487 0.9904272 0.9807487 0.9510315 0.8989985 0.81984484 0.70353085 0.52475667 0.15412536 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.4639762 0.65980405 0.7831767 0.8660702 0.9202131 0.95103025 0.9610543 0.95103025 0.9202131 0.8660702 0.7831767 0.65980405 0.4639762 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.34077275 0.57920456 0.71747255 0.80782706 0.86607176 0.8989973 0.90967524 0.8989973 0.86607176 0.80782706 0.71747255 0.57920456 0.34077275 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 6.043427f-5 0.4419233 0.612775 0.7174738 0.7831801 0.81984884 0.83167475 0.81984884 0.7831801 0.7174738 0.612775 0.4419233 6.043427f-5 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.15409686 0.44191584 0.5792005 0.6598068 0.70353055 0.71747285 0.70353055 0.6598068 0.5792005 0.44191584 0.15409686 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.34076345 0.46397528 0.5247502 0.5435158 0.5247502 0.46397528 0.34076345 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.15412536 0.2019202 0.15412536 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    simple_view(shapes[4])
    @test camera.canvas_grayscale ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.978673 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.9981043 0.9490156 0.8905761 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.9606215 0.91050285 0.8508289 0.7812944 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.9138984 0.86261004 0.80156815 0.7304569 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.90045446 0.8574364 0.8048271 0.74227124 0.66944265 0.5860398 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.83507824 0.79075146 0.7366586 0.67243195 0.5977333 0.5122531 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.75919366 0.713377 0.6576244 0.5915577 0.51482594 0.42711005 0.32812667 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.7101855 0.672366 0.6248648 0.56726426 0.49917352 0.42023253 0.33010846 0.22850767 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.61374575 0.57418245 0.52478987 0.4651388 0.39483014 0.31349003 0.22077267 0.116371654 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5057838 0.46425173 0.41274664 0.35083148 0.2780933 0.19414923 0.098636374 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.4200221 0.38592935 0.34218964 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    simple_view(shapes[5])
    @test camera.canvas_grayscale ≈ Float32[
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.81478345 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.9987888 0.0 0.84281206 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.97484577 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.94499964 0.5612857 0.8732455 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.84932184 0.0 0.8753272 0.0 0.78884614 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.06278893 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.84617716 0.094548054 0.051296055 0.0 0.0 0.0 0.5318333 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.81487197 0.0 0.7603712 0.0 0.0 0.62965536 0.0875735 0.3295365 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.38508785 0.3571475 0.0 0.0 0.0 0.0 0.0 0.0 0.36121896 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.6843042 0.0 0.43388656 0.5853236 0.44105455 0.0 0.42422578 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5340948 0.60774314 0.5798033 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
end
