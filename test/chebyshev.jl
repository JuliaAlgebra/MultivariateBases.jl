using Test
using MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Univariate" begin
    basis = maxdegree_basis(ChebyshevBasis, (x,), 3)
    @test basis.polynomials[4] == 1
    @test basis.polynomials[3] == x
    @test basis.polynomials[2] == 2x^2 - 1
    @test basis.polynomials[1] == 4x^3 - 3x
end

@testset "API degree = $degree" for degree in 0:3
    api_test(ChebyshevBasis, degree)
end
