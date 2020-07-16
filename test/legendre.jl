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

    univ_orthogonal_test(LegendreBasis, i -> 2/(2*i+1))
end

@testset "API degree = $degree" for degree in 0:3
    api_test(LegendreBasis, degree)
end

@testset "Coefficients" begin
    coefficient_test(LegendreBasis, [0.1523809523809524, 0.1523809523809524, 0.0761904761904762, -0.5714285714285714, 0.0761904761904762, -0.3428571428571428, -0.3428571428571428, 0.8] )
end
