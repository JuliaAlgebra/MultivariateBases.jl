using Test

using MultivariateBases

using DynamicPolynomials

function api_test(B::Type{<:AbstractPolynomialBasis}, degree)
    @polyvar x[1:2]
    basis = maxdegree_basis(B, x, degree)
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

@testset "Monomial" begin
    include("monomial.jl")
end
@testset "Scaled" begin
    include("scaled.jl")
end
@testset "Fixed" begin
    include("fixed.jl")
end
@testset "Chebyshev" begin
    include("chebyshev.jl")
end
