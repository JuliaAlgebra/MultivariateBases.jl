using Test
using MultivariateBases

@testset "Orthogonal" begin
    orthogonal_test(
        ProbabilistsHermiteBasis,
        x -> [1, x, x^2 - 1, x^3 - 3x, x^4 - 6x^2 + 3],
        true,
    )
    orthogonal_test(
        PhysicistsHermiteBasis,
        x -> [1, 2x, 4x^2 - 2, 8x^3 - 12x, 16x^4 - 48x^2 + 12],
        true,
    )
end

@testset "API degree = $degree" for degree in 0:3
    api_test(ProbabilistsHermiteBasis, degree)
    api_test(PhysicistsHermiteBasis, degree)
end
