using Test
using MultivariateBases
const MB = MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Linear" begin
    basis = SubBasis{MB.Monomial}([x, y])
    @test basis == SubBasis{MB.Monomial}([y, x])
    @test basis != SubBasis{MB.Monomial}([y, y^2, x])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> i^2, basis) == 4x + y
    @test coefficients(x + 4y, basis) == [4, 1]
    @test polynomial(basis[1]) == y
    @test polynomial(basis[2]) == x
    @test basis.monomials == [y, x]
    @test polynomial.(collect(basis)) == [y, x]
    @test variables(basis) == [x, y]
    @test nvariables(basis) == 2
    @test sprint(show, basis) == "SubBasis{Monomial}([y, x])"
    @test sprint(show, MIME"text/print"(), basis) ==
          "SubBasis{Monomial}([y, x])"
    @test sprint(show, MIME"text/plain"(), basis) ==
          "SubBasis{Monomial}([y, x])"
    @test sprint(print, basis) == "SubBasis{Monomial}([y, x])"
end
@testset "Affine" begin
    # It will be sorted and 1 will be moved at the end
    basis = SubBasis{MB.Monomial}([1, x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> i^2, basis) == 9x + 4y + 1
    @test coefficients(9 + x + 4y, basis) == [9, 4, 1]
end
@testset "Quadratic" begin
    basis = SubBasis{MB.Monomial}([x^2, x * y, y^2])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> i^2, basis) == 9x^2 + 4x * y + y^2
    @test coefficients(x^2 + 4x * y + 9y^2, basis) == [9, 4, 1]
    @test sprint(show, basis) == "SubBasis{Monomial}([y², xy, x²])"
    @test sprint(print, basis) == "SubBasis{Monomial}([y^2, x*y, x^2])"
end
@testset "merge_bases" begin
    basis1 = SubBasis{MB.Monomial}([x^2, y^2])
    basis2 = SubBasis{MB.Monomial}([x * y, y^2])
    @test mindegree(basis2) == 2
    @test mindegree(basis2, x) == 0
    @test mindegree(basis2, y) == 1
    @test maxdegree(basis2) == 2
    @test maxdegree(basis2, x) == 1
    @test maxdegree(basis2, y) == 2
    @test extdegree(basis2) == (2, 2)
    @test extdegree(basis2, x) == (0, 1)
    @test extdegree(basis2, y) == (1, 2)
    basis, I1, I2 = MultivariateBases.merge_bases(basis1, basis2)
    @test basis.monomials == [y^2, x * y, x^2]
    @test I1 == [1, 0, 2]
    @test I2 == [1, 2, 0]
end

@testset "API degree = $degree" for degree in 0:3
    api_test(MB.Monomial, degree)
end

@testset "Empty" begin
    basis = MB.SubBasis{MB.Monomial}(typeof(x^2)[])
    @test isempty(basis)
    @test isempty(eachindex(basis))
    p = @inferred polynomial(zeros(Int, 0, 0), basis, Int)
    @test iszero(p)
end

@testset "Enumerate" begin
    monos = [1, y, x]
    basis = MB.SubBasis{MB.Monomial}(monos)
    for (i, e) in enumerate(basis)
        @test polynomial(e) == monos[i]
    end
end

@testset "Coefficients" begin
    coefficient_test(MB.Monomial, [1, -3, 1, 1])
end
