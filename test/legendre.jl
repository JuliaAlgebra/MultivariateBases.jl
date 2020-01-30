using Test
using MultivariateBases

@testset "Orthogonal" begin
    orthogonal_test(LegendreBasis, x -> [
        1,
        x,
        (3x^2 - 1) / 2,
        (5x^3 - 3x) / 2,
        (35x^4 - 30x^2 + 3) / 8
    ], true)
end

@testset "API degree = $degree" for degree in 0:3
    api_test(LegendreBasis, degree)
end
