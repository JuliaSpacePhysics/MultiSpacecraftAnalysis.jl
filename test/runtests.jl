using MultiSpacecraftAnalysis
using Test
using Aqua

@testset "MultiSpacecraftAnalysis.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(MultiSpacecraftAnalysis)
    end
    # Write your tests here.
end
