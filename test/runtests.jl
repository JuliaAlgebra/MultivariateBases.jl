using Test

using MultivariateBases
using LinearAlgebra
using DynamicPolynomials


function api_test(B::Type{<:AbstractPolynomialBasis}, degree)
    @polyvar x[1:2]
    for basis in [
        maxdegree_basis(B, x, degree),
        basis_covering_monomials(B, monomials(x, 0:degree))
    ]
        n = binomial(2 + degree, 2)
        @test length(basis) == n
        @test typeof(copy(basis)) == typeof(basis)
        @test nvariables(basis) == 2
        @test variables(basis) == x
        @test monomialtype(typeof(basis)) == monomialtype(x[1])
        @test typeof(empty_basis(typeof(basis))) == typeof(basis)
        @test length(empty_basis(typeof(basis))) == 0
        @test polynomialtype(basis, Float64) == polynomialtype(x[1], Float64)
        @test polynomial(i -> 0.0, basis) isa polynomialtype(basis, Float64)
        @test polynomial(zeros(n, n), basis, Float64) isa polynomialtype(basis, Float64)
        @test polynomial(ones(n, n), basis, Float64) isa polynomialtype(basis, Float64)
    end
end

function univ_orthogonal_test(B::Type{<:AbstractMultipleOrthogonalBasis}, univ::Function; kwargs...)
    @polyvar x
    basis = maxdegree_basis(B, [x], 4)
    for i = 1:length(basis)
        @test isapprox(dot(basis[i], basis[i], B), univ(maxdegree(basis[i])); kwargs...)
        for j = 1:i-1
            @test isapprox(dot(basis[i], basis[j], B), 0.0; kwargs...)
        end
    end
end

function orthogonal_test(B::Type{<:AbstractMultipleOrthogonalBasis}, univ::Function, even_odd_separated::Bool)
    @test MultivariateBases.even_odd_separated(B) == even_odd_separated
    @polyvar x y
    univariate_x = univ(x)
    univariate_y = univ(y)

    @testset "Univariate $var" for (var, univ) in [(x, univariate_x), (y, univariate_y)]
        basis = maxdegree_basis(B, (var,), 4)
        for (i, e) in enumerate(basis)
            @test e == univ[(length(univ) +1 - i)]
        end
    end

    @testset "basis_covering_monomials" begin
        monos = basis_covering_monomials(B, monovec([x^2 * y, y^2]))
        if even_odd_separated
            exps = [(2, 1), (0, 2), (0, 1), (0, 0)]
        else
            exps = [(2, 1), (2, 0), (1, 1), (0, 2), (1, 0), (0, 1), (0, 0)]
        end
        for (i, e) in enumerate(monos)
            @test e == univariate_x[exps[i][1] + 1] * univariate_y[exps[i][2] + 1]
        end
        monos = basis_covering_monomials(B, monovec([x^4, x^2, x]))
        if even_odd_separated
            exps = [4, 2, 1, 0]
        else
            exps = [4, 3, 2, 1, 0]
        end
        for (i, e) in enumerate(monos)
            @test e == univariate_x[exps[i] + 1]
        end
    end
end

function coefficient_test(B::Type{<:AbstractPolynomialBasis}, coefs; kwargs...)
    @polyvar x y
    p = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1
    cc = coefficients(p, B)
    for (i, c) in enumerate(cc)
        @test isapprox(c, coefs[i]; kwargs...)
    end

    mons = basis_covering_monomials(B, monomials(p))
    @test isapprox(p, polynomial(cc, mons); kwargs...)
    
end

@testset "Monomial" begin
    include("monomial.jl")
end
@testset "Scaled" begin
    include("scaled.jl")
end
@testset "Fixed" begin
    include("fixed.jl")
end
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
