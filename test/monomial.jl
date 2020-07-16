using Test
using MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Linear" begin
    basis = MonomialBasis([x, y])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(i -> i^2, basis) == x + 4y
    @test coefficients(x + 4y, MonomialBasis) == [1, 4]
end
@testset "Affine" begin
    # It will be sorted and 1 will be moved at the end
    basis = MonomialBasis([1, x, y])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(i -> i^2, basis) == x + 4y + 9
    @test coefficients(9 + x + 4y, MonomialBasis) == [1, 4, 9]
end
@testset "Quadratic" begin
    basis = MonomialBasis([x^2, x*y, y^2])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(i -> i^2, basis) == x^2 + 4x*y + 9y^2
    @test coefficients(x^2 + 4x*y + 9y^2, MonomialBasis) == [1, 4, 9]
end
@testset "API degree = $degree" for degree in 0:3
    api_test(MonomialBasis, degree)
end

@testset "Enumerate" begin
    monos = [1, y, x]
    basis = MonomialBasis(monos)
    for (i, e) in enumerate(basis)
        @test e == monos[length(monos)+1-i]
    end
end

@testset "Coefficients" begin
    coefficient_test(MonomialBasis, [1, 1, -3, 1])
end
