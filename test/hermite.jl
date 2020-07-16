using Test
using MultivariateBases

@testset "Orthogonal" begin
    orthogonal_test(ProbabilistsHermiteBasis, x -> [
        1,
        x,
        x^2 - 1,
        x^3 - 3x,
        x^4 - 6x^2 + 3
    ], true)
    univ_orthogonal_test(ProbabilistsHermiteBasis, i-> √(2*π)*factorial(i))
    orthogonal_test(PhysicistsHermiteBasis, x -> [
        1,
        2x,
        4x^2 - 2,
        8x^3 - 12x,
        16x^4 - 48x^2 + 12
    ], true)
    univ_orthogonal_test(PhysicistsHermiteBasis, i-> √(π)*factorial(i)*2^i, atol = 1e-12) #precision issue
end

@testset "API degree = $degree" for degree in 0:3
    api_test(ProbabilistsHermiteBasis, degree)
    api_test(PhysicistsHermiteBasis, degree)
end

@testset "Coefficients" begin
    coefficient_test(ProbabilistsHermiteBasis, [1.0, 1.0, 1.0, 9.0, 1.0, 6.0, 6.0, 4.0]; atol = 1e-12)
    coefficient_test(PhysicistsHermiteBasis, [0.015625, 0.015625, 0.03125, 0.1875, 0.03125, 0.1875, 0.1875, 1.0])
end
