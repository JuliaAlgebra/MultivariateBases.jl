using Test

using MultivariateBases

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
        @test monomial_type(typeof(basis)) == monomial_type(x[1])
        @test typeof(empty_basis(typeof(basis))) == typeof(basis)
        @test length(empty_basis(typeof(basis))) == 0
        @test polynomial_type(basis, Float64) == polynomial_type(x[1], Float64)
        @test polynomial(i -> 0.0, basis) isa polynomial_type(basis, Float64)
        @test polynomial(zeros(n, n), basis, Float64) isa polynomial_type(basis, Float64)
        @test polynomial(ones(n, n), basis, Float64) isa polynomial_type(basis, Float64)
    end
end

function orthogonal_test(B::Type{<:AbstractMultipleOrthogonalBasis}, univ::Function, even_odd_separated::Bool)
    @test MultivariateBases.even_odd_separated(B) == even_odd_separated
    @polyvar x y
    univariate_x = univ(x)
    univariate_y = univ(y)

    @testset "Univariate $var" for (var, univ) in [(x, univariate_x), (y, univariate_y)]
        basis = maxdegree_basis(B, (var,), 4)
        for i in 1:5
            @test basis.polynomials[length(basis) + 1 - i] == univ[i]
        end
    end

    @testset "basis_covering_monomials" begin
        monos = basis_covering_monomials(B, monomial_vector([x^2 * y, y^2]))
        if even_odd_separated
            exps = [(2, 1), (0, 2), (0, 1), (0, 0)]
        else
            exps = [(2, 1), (2, 0), (1, 1), (0, 2), (1, 0), (0, 1), (0, 0)]
        end
        for i in 1:length(monos)
            @test monos.polynomials[i] == univariate_x[exps[i][1] + 1] * univariate_y[exps[i][2] + 1]
        end
        monos = basis_covering_monomials(B, monomial_vector([x^4, x^2, x]))
        if even_odd_separated
            exps = [4, 2, 1, 0]
        else
            exps = [4, 3, 2, 1, 0]
        end
        for i in 1:length(monos)
            @test monos.polynomials[i] == univariate_x[exps[i] + 1]
        end
    end
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
