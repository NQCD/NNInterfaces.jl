using NNInterfaces
using Test
using FiniteDiff

function finite_difference_gradient(model, R)
    f(x) = potential(model, x)
    FiniteDiff.finite_difference_gradient(f, R)
end

@testset "NNInterfaces.jl" begin

    @testset "H2AgModel" begin
        model = H2AgModel()

        R = [
            0.313365  0.783967;
            0.527937  0.106836;
            0.909677  0.211681
        ]

        Vtest = 0.43294117958582595

        Dtest = [
            0.193842   -0.154042;
            -0.144525    0.169208;
            -0.0823145   0.31384
        ]

        Ftest = [
            0.0877215    0.0110317    0.0425049    0.0342595  -0.00260693   0.0175242;
            0.0110317    0.140706     0.00210935  -0.0432715  -0.0104426    0.0534511;
            0.0425049    0.00210935   0.0159194   -0.0638902  -0.0364848   -0.013972;
            0.0342595   -0.0432715   -0.0638902    0.0648546   0.0350274   -0.0230929;
            -0.00260693  -0.0104426   -0.0364848    0.0350274   0.0379259    0.0194911;
            0.0175242    0.0534511   -0.013972    -0.0230929   0.0194911    0.0461783
        ]

        @test potential(model, R) ≈ Vtest rtol=1e-3
        @test derivative(model, R) ≈ Dtest rtol=1e-3
        @test finite_difference_gradient(model, R) ≈ Dtest rtol=1e-3
        @test friction(model, R) ≈ Ftest rtol=1e-3
    end

end
