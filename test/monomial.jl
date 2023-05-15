using Test
using MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Linear" begin
    basis = MonomialBasis([x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> i^2, basis) == 4x + y
    @test coefficients(x + 4y, MonomialBasis) == [4, 1]
    @test basis[1] == y
    @test basis[2] == x
end
@testset "Affine" begin
    # It will be sorted and 1 will be moved at the end
    basis = MonomialBasis([1, x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> i^2, basis) == 9x + 4y + 1
    @test coefficients(9 + x + 4y, MonomialBasis) == [9, 4, 1]
end
@testset "Quadratic" begin
    basis = MonomialBasis([x^2, x * y, y^2])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> i^2, basis) == x^2 + 4x * y + 9y^2
    @test coefficients(x^2 + 4x * y + 9y^2, MonomialBasis) == [9, 4, 1]
end
@testset "API degree = $degree" for degree in 0:3
    api_test(MonomialBasis, degree)
end
