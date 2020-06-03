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
    @test generators(basis) == [y, x]
    @test collect(basis) == [y, x]
    @test length(basis) == 2
    @test firstindex(basis) == 1
    @test lastindex(basis) == 2
    @test variables(basis) == [x, y]
    @test nvariables(basis) == 2
    @test sprint(show, basis) == "MonomialBasis([y, x])"
    @test sprint(show, MIME"text/print"(), basis) == "MonomialBasis([y, x])"
    @test sprint(show, MIME"text/plain"(), basis) == "MonomialBasis([y, x])"
    @test sprint(print, basis) == "MonomialBasis([y, x])"
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
    @test polynomial(i -> i^2, basis) == 9x^2 + 4x * y + y^2
    @test coefficients(x^2 + 4x * y + 9y^2, MonomialBasis) == [9, 4, 1]
    @test sprint(show, basis) == "MonomialBasis([y², xy, x²])"
    @test sprint(print, basis) == "MonomialBasis([y^2, x*y, x^2])"
end
@testset "merge_bases" begin
    basis1 = MonomialBasis([x^2, y^2])
    basis2 = MonomialBasis([x * y, y^2])
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
    api_test(MonomialBasis, degree)
end

@testset "Empty" begin
    basis = MonomialBasis(typeof(x^2)[])
    @test isempty(basis)
    @test isempty(eachindex(basis))
    p = @inferred polynomial(zeros(Int, 0, 0), basis, Int)
    @test iszero(p)
end

@testset "Enumerate" begin
    monos = [1, y, x]
    basis = MonomialBasis(monos)
    for (i, e) in enumerate(basis)
        @test e == monos[i]
    end
end

@testset "Coefficients" begin
    coefficient_test(MonomialBasis, [1, -3, 1, 1])
end
