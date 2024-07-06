using Test

import MutableArithmetics as MA
import StarAlgebras as SA
using MultivariateBases
const MB = MultivariateBases
using LinearAlgebra
using DynamicPolynomials

function _test_op(op, args...)
    result = @inferred op(args...)
    @test typeof(result) == MA.promote_operation(op, typeof.(args)...)
    return result
end

function _test_basis(basis)
    B = typeof(basis)
    @test typeof(MB.algebra(basis)) == MA.promote_operation(MB.algebra, B)
    @test typeof(MB.constant_algebra_element(B, 1)) ==
          MB.constant_algebra_element_type(B, Int)
end

struct TypeA end
struct TypeB end
Base.zero(::Type{TypeA}) = TypeA()
Base.:*(::Float64, ::TypeA) = TypeB()

function api_test(B::Type{<:MB.AbstractMonomialIndexed}, degree)
    @polyvar x[1:2]
    M = typeof(prod(x))
    full_basis = FullBasis{B,M}()
    _test_basis(full_basis)
    @test sprint(show, MB.algebra(full_basis)) ==
          "Polynomial algebra of $B basis"
    for basis in [
        maxdegree_basis(full_basis, x, degree),
        explicit_basis_covering(
            full_basis,
            MB.SubBasis{MB.Monomial}(monomials(x, 0:degree)),
        ),
        explicit_basis_covering(
            full_basis,
            MB.SubBasis{ScaledMonomial}(monomials(x, 0:degree)),
        ),
    ]
        _test_basis(basis)
        @test basis isa MB.explicit_basis_type(typeof(full_basis))
        for i in eachindex(basis)
            mono = basis.monomials[i]
            poly = MB.Polynomial{B}(mono)
            @test basis[i] == poly
            @test basis[poly] == i
        end
        n = binomial(2 + degree, 2)
        @test length(basis) == n
        @test firstindex(basis) == 1
        @test lastindex(basis) == n
        @test typeof(copy(basis)) == typeof(basis)
        @test nvariables(basis) == 2
        @test variables(basis) == x
        @test monomial_type(typeof(basis)) == monomial_type(x[1])
        @test typeof(empty_basis(typeof(basis))) == typeof(basis)
        @test length(empty_basis(typeof(basis))) == 0
        @test polynomial_type(basis, Float64) == polynomial_type(x[1], Float64)
        #@test polynomial(i -> 0.0, basis) isa polynomial_type(basis, Float64)
        a = MB.algebra_element(ones(length(basis)), basis)
        _test_op(MB.implicit, a)
        @test SA.star(a) == a
    end
    mono = x[1]^2 * x[2]^3
    p = MB.Polynomial{B}(mono)
    @test full_basis[p] == mono
    @test full_basis[mono] == p
    @test polynomial_type(mono, B == Monomial ? TypeA : TypeB) ==
          polynomial_type(typeof(p), TypeA)
    a = MB.algebra_element(p)
    @test variables(a) == x
    @test typeof(polynomial(a)) == polynomial_type(typeof(a))
    @test typeof(polynomial(a)) == polynomial_type(typeof(p), Int)
    @test a ≈ a
    if B == MB.Monomial
        @test a ≈ p.monomial
        @test p.monomial ≈ a
    else
        @test !(a ≈ p.monomial)
        @test !(p.monomial ≈ a)
    end
    _wrap(s) = (B == MB.Monomial ? s : "$B($s)")
    @test sprint(show, p) == _wrap(sprint(show, p.monomial))
    @test sprint(print, p) == _wrap(sprint(print, p.monomial))
    mime = MIME"text/latex"()
    @test sprint(show, mime, p) ==
          "\$\$ " *
          _wrap(MB.SA.trim_LaTeX(mime, sprint(show, mime, p.monomial))) *
          " \$\$"
    const_mono = constant_monomial(prod(x))
    const_poly = MB.Polynomial{B}(const_mono)
    const_alg_el = MB.algebra_element(const_poly)
    for other in (const_mono, 1, const_alg_el)
        @test _test_op(+, other, const_alg_el) ≈ _test_op(*, 2, other)
        @test _test_op(+, const_alg_el, other) ≈ _test_op(*, 2, other)
        @test iszero(_test_op(-, other, const_alg_el))
        @test iszero(_test_op(-, const_alg_el, other))
    end
    @test iszero(MA.operate!(-, polynomial(const_alg_el), const_alg_el))
    @test iszero(MA.operate!(+, -polynomial(const_alg_el), const_alg_el))
    @test typeof(MB.sparse_coefficients(sum(x))) ==
          MA.promote_operation(MB.sparse_coefficients, typeof(sum(x)))
    @test typeof(MB.algebra_element(sum(x))) ==
          MA.promote_operation(MB.algebra_element, typeof(sum(x)))
    @test const_alg_el(x => ones(length(x))) == const_mono
    @test const_alg_el(ones(length(x))) == const_mono
    @test const_alg_el(ones(length(x))...) == const_mono
    @test const_alg_el(tuple(ones(length(x))...)) == const_mono
    @test subs(const_alg_el, x => ones(length(x))) == const_mono
    @test const_alg_el(x[1] => x[2], x[2] => x[1]) == const_mono
    @test differentiate(const_alg_el, x) == differentiate(const_mono, x)
    @test differentiate(const_alg_el, x, 2) == differentiate(const_mono, x, 2)
    @test MA.operate_to!(polynomial(const_mono), *, const_alg_el, const_mono) ==
          const_mono
    @test MA.operate_to!(polynomial(const_mono), *, const_mono, const_alg_el) ==
          const_mono
    @test dot([const_alg_el], [const_mono]) == const_mono
    @test dot([const_mono], [const_alg_el]) == const_mono
end

function univ_orthogonal_test(
    B::Type{<:AbstractMultipleOrthogonal},
    univ::Function;
    kwargs...,
)
    @polyvar x
    basis = maxdegree_basis(FullBasis{B,monomial_type(x)}(), [x], 4)
    for i in eachindex(basis)
        p_i = polynomial(basis[i])
        @test isapprox(dot(p_i, p_i, B), univ(maxdegree(p_i)); kwargs...)
        for j in 1:i-1
            @test isapprox(dot(p_i, polynomial(basis[j]), B), 0.0; kwargs...)
        end
    end
end

function orthogonal_test(
    B::Type{<:AbstractMultipleOrthogonal},
    univ::Function,
    even_odd_separated::Bool,
)
    @test MultivariateBases.even_odd_separated(B) == even_odd_separated
    @polyvar x y
    univariate_x = univ(x)
    univariate_y = univ(y)

    @testset "Univariate $var" for (var, univ) in
                                   [(x, univariate_x), (y, univariate_y)]
        basis = maxdegree_basis(FullBasis{B,monomial_type(var)}(), (var,), 4)
        for i in 1:5
            @test polynomial(basis[i]) == univ[i]
        end
    end

    @testset "explicit_basis_covering" begin
        basis = explicit_basis_covering(
            FullBasis{B,typeof(x * y)}(),
            SubBasis{MB.Monomial}(monomial_vector([x^2 * y, y^2])),
        )
        if even_odd_separated
            exps = [(0, 0), (0, 1), (0, 2), (2, 1)]
        else
            exps = [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0), (2, 1)]
        end
        for i in eachindex(basis)
            @test polynomial(basis[i]) ==
                  univariate_x[exps[i][1]+1] * univariate_y[exps[i][2]+1]
        end
        basis = explicit_basis_covering(
            FullBasis{B,typeof(x^2)}(),
            SubBasis{MB.Monomial}(monomial_vector([x^4, x^2, x])),
        )
        if even_odd_separated
            exps = [0, 1, 2, 4]
        else
            exps = 0:4
        end
        for i in eachindex(basis)
            @test polynomial(basis[i]) == univariate_x[exps[i]+1]
        end
    end
end

function coefficient_test(basis::SubBasis, p, coefs; kwargs...)
    cc = coefficients(p, basis)
    @test isapprox(coefs, cc; kwargs...)
    @test isapprox(p, algebra_element(cc, basis); kwargs...)
end

function coefficient_test(basis::SubBasis, coefs::AbstractVector; kwargs...)
    return coefficient_test(
        basis,
        sum(generators(basis) .* coefs),
        coefs;
        kwargs...,
    )
end

function coefficient_test(
    B::Type{<:MB.AbstractMonomialIndexed},
    coefs;
    kwargs...,
)
    @polyvar x y
    p = x^4 * y^2 + x^2 * y^4 - 3 * x^2 * y^2 + 1
    basis = explicit_basis_covering(
        FullBasis{B,typeof(x * y)}(),
        SubBasis{MB.Monomial}(monomials(p)),
    )
    coefficient_test(basis, p, coefs; kwargs...)
    return
end

@testset "Monomial" begin
    include("monomial.jl")
end
@testset "Scaled" begin
    include("scaled.jl")
end
#@testset "Fixed" begin
#    include("fixed.jl")
#end
#@testset "Orthonormal" begin
#    include("orthonormal.jl")
#end
@testset "Hermite" begin
    include("hermite.jl")
end
@testset "Laguerre" begin
    include("laguerre.jl")
end
@testset "Legendre" begin
    include("legendre.jl")
end
@testset "Chebyshev" begin
    include("chebyshev.jl")
end
@testset "Quotient" begin
    include("quotient.jl")
end
