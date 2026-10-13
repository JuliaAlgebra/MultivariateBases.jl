module TestSharedCore

using Test
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import MultivariateBases as MB
import StarAlgebras as SA
import LinearAlgebra

function indexed_element(basis, exp, coeff)
    return MB.algebra_element(
        SA.SparseCoefficients([exp], [coeff], SA.comparable(basis)),
        basis,
    )
end

DP.@polyvar x y

@testset "Shared monomial core" begin
    @test MB.Polynomial === MP.Polynomial
    @test MB.FullBasis === MP.FullBasis
    @test MB.MStruct === MP.MStruct
    full = MB.FullBasis{MB.Monomial}(MP.variables(x * y))
    a = indexed_element(full, [1, 0], 2)
    b = indexed_element(full, [0, 1], 3)
    @test collect(SA.values(SA.coeffs(a * b))) == [6]
    @test only(SA.keys(SA.coeffs(a * b))) == [1, 1]
    @test MB.algebra_element(a) === a
    @test MB.Polynomial{MB.Monomial}(x) isa MP.AbstractMonomial
end

@testset "Scaled basis" begin
    full = MB.FullBasis{MB.ScaledMonomial}(MP.variables(x * y))
    @test MB.Monomial !== MB.ScaledMonomial
    @test one(full.map) isa MB.Polynomial{MB.ScaledMonomial}
    a = indexed_element(full, [1, 0], 1.0)
    b = indexed_element(full, [0, 1], 1.0)
    @test only(SA.values(SA.coeffs(a * b))) ≈ inv(sqrt(2.0))
    mono = MB.FullBasis{MB.Monomial}(MP.variables(x * y))
    expanded = MB.convert_basis(mono, a * b)
    @test only(SA.values(SA.coeffs(expanded))) ≈ 1.0
    @test only(SA.keys(SA.coeffs(expanded))) == [1, 1]
    @test MB.algebra_element(full[[1, 1]]) isa MP.AbstractPolynomial
end

@testset "Chebyshev basis" begin
    full = MB.FullBasis{MB.Chebyshev}(MP.variables(x))
    a = indexed_element(full, [1], 1 // 1)
    squared = a * a
    @test collect(SA.keys(SA.coeffs(squared))) == [[0], [2]]
    @test collect(SA.values(SA.coeffs(squared))) == [1 // 2, 1 // 2]
    mono = MB.FullBasis{MB.Monomial}(MP.variables(x))
    expanded = MB.convert_basis(mono, squared)
    @test collect(SA.keys(SA.coeffs(expanded))) == [[2]]
    @test collect(SA.values(SA.coeffs(expanded))) == [1 // 1]
    @test squared(x => 2) == 4
    @test MP.differentiate(squared, x) == 2 * x
    t2 = indexed_element(full, [2], 1 // 1)
    @test LinearAlgebra.dot(t2, t2, MB.Chebyshev) ≈ π / 2
    @test LinearAlgebra.dot(t2, one(t2), MB.Chebyshev) ≈ 0 atol = 1e-12
    sub = SA.SubBasis(full, [[0], [1], [2]])
    @test MP.coefficients(t2, sub) ≈ [0, 0, 1]
    @test MP.coefficients(x^2, sub) ≈ [1 / 2, 0, 1 / 2]
end

end # module
