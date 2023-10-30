import Rays

@testset "blur_kernel" begin
    kernel, Δi_max = Rays.get_blur_kernel(1.6)
    @test sum(kernel) ≈ 1.0
    @test all(kernel .>= 0.0)
    @test kernel ≈ Float32[0.00029104948, 0.2252549, 0.5489081, 0.2252549, 0.00029104948]
end
