using Test
using MultivariateBases
const MB = MultivariateBases
using MultivariatePolynomials

function test_monomial(x, y)
    @testset "StarAlgebras" begin
        vars = MB.Variables{MB.Monomial}(variables(x))
        @test vars(exponents(x^2)) == vars(exponents(x^2))
        @test vars(exponents(x^3)) != vars(exponents(x^2))
        o = vars(exponents(constant_monomial(x^2)))
        @test isone(o)
    end

    @testset "Linear" begin
        basis = SubBasis{MB.Monomial}([x, y])
        @test basis == SubBasis{MB.Monomial}([y, x])
        @test basis != SubBasis{MB.Monomial}([y, y^2, x])
        @test polynomial_type(basis, Int) == polynomial_type(x * y, Int)
        @test polynomial(i -> i^2, basis) == 4x + y
        @test coefficients(x + 4y, basis) == [4, 1]
        @test polynomial(basis[1]) == y
        @test polynomial(basis[2]) == x
        @test MB.keys_as_monomials(basis) == [y, x]
        @test polynomial.(collect(basis)) == [y, x]
        @test variables(basis) == variables(x * y)
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
        @test polynomial_type(basis, Int) == polynomial_type(x * y, Int)
        poly = 9x + 4y + 1
        @test polynomial(i -> i^2, basis) == poly
        @test coefficients(9 + x + 4y, basis) == [9, 4, 1]
        @test MB.explicit_basis(poly) == basis
    end

    @testset "Quadratic" begin
        basis = SubBasis{MB.Monomial}([x^2, x * y, y^2])
        @test polynomial_type(basis, Int) == polynomial_type(x * y, Int)
        @test polynomial(i -> i^2, basis) == 9x^2 + 4x * y + y^2
        @test coefficients(x^2 + 4x * y + 9y^2, basis) == [9, 4, 1]
        @test sprint(show, basis) == "SubBasis{Monomial}([y², xy, x²])"
        @test sprint(print, basis) == "SubBasis{Monomial}([y^2, x*y, x^2])"
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

    @testset "promote_bases_with_maps" begin
        a = MB.algebra_element(x - x^2)
        b = MB.FullBasis{MB.Monomial}(x * y)
        a2, b2 = SA.promote_bases(a, b)
        @test collect.(a2.coeffs.basis_elements) == [[1, 0], [2, 0]]
        @test collect.(a2.coeffs.basis_elements) == [[1, 0], [2, 0]]
        @test SA.basis(a) != b
        @test SA.basis(a2) == b
        b2 === b
    end

    @testset "hash" begin
        monosx = [1, x]
        basisx = MB.SubBasis{MB.Monomial}(monosx)
        scaled_basisx = MB.SubBasis{MB.ScaledMonomial}(monosx)
        monosy = [1, y]
        basisy = MB.SubBasis{MB.Monomial}(monosy)
        @test hash(basisx[2]) != hash(scaled_basisx[2])
        @test hash(basisx[2]) != hash(basisy[2])
    end
end

import DynamicPolynomials
test_monomial((DynamicPolynomials.@polyvar x y)...)

import TypedPolynomials
test_monomial((TypedPolynomials.@polyvar x y)...)

import MutableArithmetics as MA

DynamicPolynomials.@ncpolyvar x y
p = x * y + x
q = y + x * y
a = MB.algebra_element(p)
b = MB.algebra_element(q)
c = MB.algebra_element(p * q)
MA.operate!(zero, c)
MA.operate_to!(c, *, a, b)
@edit a * b
p * q
