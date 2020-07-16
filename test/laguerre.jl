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

    univ_orthogonal_test(LaguerreBasis, i -> 1; atol = 1e-12) # there are precision issues
end

@testset "API degree = $degree" for degree in 0:3
    api_test(LaguerreBasis, degree)
end

@testset "Coefficients" begin
    coefficient_test(LaguerreBasis, [48.0, 48.0, -96.0, -192.0, -192.0, -96.0, 48.0, 384.0, 564.0, 384.0, 48.0, -192.0, -744.0, -744.0, -192.0, 324.0, 720.0, 324.0, -264.0, -264.0, 85.0])
end
