using Test
using MultivariateBases

@testset "Orthogonal" begin
    orthogonal_test(ChebyshevBasis, x -> [
        1,
        x,
        2x^2 - 1,
        4x^3 - 3x,
        8x^4 - 8x^2 + 1
    ], true)
    univ_orthogonal_test(ChebyshevBasis, i -> π*(1/2 + (0^i)/2), atol = 1e-12)
    orthogonal_test(ChebyshevBasisSecondKind, x -> [
        1,
        2x,
        4x^2 - 1,
        8x^3 - 4x,
        16x^4 - 12x^2 + 1
    ], true)

    univ_orthogonal_test(ChebyshevBasisSecondKind, i -> π/2, atol = 1e-12)
end

@testset "API degree = $degree" for degree in 0:3
    api_test(ChebyshevBasis, degree)
    api_test(ChebyshevBasisSecondKind, degree)
end


@testset "Coefficients" begin
    coefficient_test(ChebyshevBasis,[0.0625, 0.0625, 0.0625, -0.25, 0.0625, -0.3125, -0.3125, 0.625])
    coefficient_test(ChebyshevBasisSecondKind, [0.015625, 0.015625, 0.015625, -0.09375, 0.015625, -0.109375, -0.109375, 0.875])
end
