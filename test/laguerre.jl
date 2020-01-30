using Test
using MultivariateBases

@testset "Orthogonal" begin
    orthogonal_test(LaguerreBasis, x -> [
        1,
        1 - x,
        (x^2 - 4x + 2) / 2,
        (-x^3 + 9x^2 - 18x + 6) / 6,
        (x^4 - 16x^3 + 72x^2 - 96x + 24) / 24
    ], false)
end

@testset "API degree = $degree" for degree in 0:3
    api_test(LaguerreBasis, degree)
end
