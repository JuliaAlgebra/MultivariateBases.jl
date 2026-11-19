using Test
using DynamicPolynomials
import MultivariateBases as MB
import StarAlgebras as SA

import SparseArrays

function _test_quotient(matrix)
    @polyvar x y
    sub = MB.SubBasis{MB.Monomial}([1, y, x])
    @test_throws DimensionMismatch MB.QuotientBasis(sub, ones(2, 2))
    # Quotient by the span of `x - y`
    basis = MB.QuotientBasis(sub, matrix)
    @test length(basis) == 2
    @test sprint(show, basis) == "Quotient of dimension 2 of the span of $sub"
    @test coefficients(2 + 3x + y, basis) == [2, 4]
    @test iszero(coefficients(x - y, basis))
    @test_throws ArgumentError coefficients(x^2, basis)
    a = MB.algebra_element(2 + 3x + y)
    @test SA.coeffs(a, basis) == [2, 4]
    explicit = MB.explicit_basis(a)
    exp = MB.algebra_element(SA.coeffs(a, explicit), explicit)
    @test SA.coeffs(exp, basis) == [2, 4]
    @test SA.adjoint_coeffs(
        [5, 7],
        MB.SubBasis{MB.Monomial}([y, x^2]),
        basis,
    ) == [7, 0]
    q = MB.algebra_element([5, 7], basis)
    @test SA.coeffs(q) == [5, 7]
    return
end

@testset "QuotientBasis" begin
    _test_quotient([1 0 0; 0 1 1])
    _test_quotient(SparseArrays.sparse([1 0 0; 0 1 1]))
end

@testset "QuotientBasis modulo x^2 + 1" begin
    @polyvar x y
    sub = MB.SubBasis{MB.Monomial}([1, y, x, x^2, x^2 * y, x^3])
    # Remainders `-1`, `-y` and `-x` of `x^2`, `x^2 * y` and `x^3`
    R = [
        1 0 0 -1 0 0
        0 1 0 0 -1 0
        0 0 1 0 0 -1
    ]
    basis = MB.QuotientBasis(sub, SparseArrays.sparse(R))
    @test length(basis) == 3
    p = x^3 - 2x^2 * y + 3x^2
    @test MB.keys_as_monomials(sub) == [1, y, x, x^2, x^2 * y, x^3]
    @test coefficients(p, basis) == [-3, 2, -1]
    @test SA.adjoint_coeffs([1, 2, 3], sub, basis) == [1, 2, 3, -1, -2, -3]
end
