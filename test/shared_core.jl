module TestSharedCore

using Test
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import MultivariateBases as MB
import StarAlgebras as SA
import MutableArithmetics as MA
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

    @testset "Scaled term products stay terms" begin
        alg = parent(a)
        t, u = SA.Term(alg, [1, 0], 2), SA.Term(alg, [0, 1], 3)
        product = @inferred t * u
        @test product isa SA.Term{Float64}
        @test product.index == [1, 1]
        @test SA.coefficient(product) ≈ 6 / sqrt(2)
        @test (@inferred MA.promote_operation(*, typeof(t), typeof(u))) ===
              typeof(product)
        z = @inferred zero(t) * u
        @test z isa typeof(product)
        @test iszero(z)
        @test (@inferred t * b) isa SA.AlgebraElement{Float64}
        @test (@inferred a * u) isa SA.AlgebraElement{Float64}
        for op in (MA.add_mul, MA.sub_mul), (left, right) in ((t, b), (a, u))
            f = a + b
            expected = op(f, left, right)
            @test (@inferred MA.operate!!(op, f, left, right)) === f
            @test f == expected
        end

        # Promote disjoint sets of variables before multiplying the basis keys.
        bx = MB.FullBasis{MB.ScaledMonomial}(MP.variables(x))
        by = MB.FullBasis{MB.ScaledMonomial}(MP.variables(y))
        tx = SA.Term(parent(indexed_element(bx, [1], 2)), [1], 2)
        ty = SA.Term(parent(indexed_element(by, [1], 3)), [1], 3)
        promoted = @inferred tx * ty
        @test promoted isa typeof(product)
        @test SA.basis(parent(promoted)) == full
        @test promoted.index == product.index
        @test SA.coefficient(promoted) ≈ SA.coefficient(product)
        @test (@inferred MA.promote_operation(*, typeof(tx), typeof(ty))) ===
              typeof(promoted)
        product.index[1] = 5
        @test t.index == [1, 0]
        @test u.index == [0, 1]
        @test tx.index == ty.index == [1]
    end
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

    @testset "General term products are algebra elements" begin
        alg = parent(a)
        for i in (0, 1), c in (0, 3)
            t, u = SA.Term(alg, [i], c), SA.Term(alg, [1], 2)
            product = @inferred t * u
            @test product isa SA.AlgebraElement{Rational{Int}}
            @test product == SA.algebra_element(t) * SA.algebra_element(u)
            @test (@inferred MA.promote_operation(*, typeof(t), typeof(u))) ===
                  typeof(product)
        end
    end
end

end # module
