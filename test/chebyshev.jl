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
    orthogonal_test(ChebyshevBasisSecondKind, x -> [
        1,
        2x,
        4x^2 - 1,
        8x^3 - 4x,
        16x^4 - 12x^2 + 1
    ], true)
end

@testset "API degree = $degree" for degree in 0:3
    api_test(ChebyshevBasis, degree)
    api_test(ChebyshevBasisSecondKind, degree)
end
