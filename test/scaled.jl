using Test
using MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Linear" begin
    basis = ScaledMonomialBasis([x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Float64)
    @test polynomial(i -> i^2, basis) == x + 4y
    @test coefficients(x + 4y, ScaledMonomialBasis) == [1, 4]
end
@testset "Affine" begin
    # It will be sorted and 1 will be moved at the end
    basis = ScaledMonomialBasis([1, x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Float64)
    @test polynomial(i -> i^2, basis) == x + 4y + 9
    @test coefficients(9 + x + 4y, ScaledMonomialBasis) == [1, 4, 9]
end
@testset "Quadratic" begin
    basis = ScaledMonomialBasis([x^2, x * y, y^2])
    @test polynomial_type(basis, Int) == polynomial_type(x, Float64)
    @test polynomial(i -> i^2, basis) == x^2 + 4 * √2 * x * y + 9y^2
    @test coefficients(x^2 + 4x * y + 9y^2, ScaledMonomialBasis) ==
          [1, 4 / √2, 9]
end
@testset "API degree = $degree" for degree in 0:3
    api_test(ScaledMonomialBasis, degree)
end
