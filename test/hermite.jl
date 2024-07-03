using Test
using MultivariateBases

@testset "Orthogonal" begin
    orthogonal_test(
        MB.ProbabilistsHermite,
        x -> [1, x, x^2 - 1, x^3 - 3x, x^4 - 6x^2 + 3],
        true,
    )
    univ_orthogonal_test(MB.ProbabilistsHermite, i -> √(2 * π) * factorial(i))
    orthogonal_test(
        MB.PhysicistsHermite,
        x -> [1, 2x, 4x^2 - 2, 8x^3 - 12x, 16x^4 - 48x^2 + 12],
        true,
    )
    univ_orthogonal_test(
        MB.PhysicistsHermite,
        i -> √(π) * factorial(i) * 2^i;
        atol = 1e-12,
    ) #precision issue
end

@testset "API degree = $degree" for degree in 0:3
    api_test(MB.ProbabilistsHermite, degree)
    api_test(MB.PhysicistsHermite, degree)
end

@testset "Coefficients" begin
    @polyvar x
    coefficient_test(
        MB.ProbabilistsHermite,
        [4, 6, 6, 1, 9, 1, 1, 1];
        atol = 1e-12,
    )
    coefficient_test(
        MB.PhysicistsHermite,
        reverse([
            0.015625,
            0.015625,
            0.03125,
            0.1875,
            0.03125,
            0.1875,
            0.1875,
            1.0,
        ]),
    )
    M = typeof(x^2)
    mono = MB.FullBasis{MB.Monomial,M}()
    basis = MB.FullBasis{MB.PhysicistsHermite,M}()
    err = ErrorException("Convertion from `$mono` to `$basis` not implemented yet")
    @test_throws err SA.coeffs(MB.algebra_element(x + 1), basis)
end
