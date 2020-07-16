using Test
using MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Polynomials" begin
    basis = FixedPolynomialBasis([1 - x^2, x^2 + 2])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(one, basis) == 3
end
@testset "Monomial" begin
    basis = FixedPolynomialBasis([x, x^2])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(i -> i^2, basis) == 4x^2 + x
end
@testset "Terms" begin
    basis = FixedPolynomialBasis([1, x, x^2])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(i -> (-1)^i, basis) == -x^2 + x - 1
end
@testset "Linear" begin
    basis = FixedPolynomialBasis([x, y])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(identity, basis) == x + 2y
end
@testset "One variable" begin
    basis = FixedPolynomialBasis([x])
    @test polynomialtype(basis, Int) == polynomialtype(x, Int)
    @test polynomial(i -> 5, basis) == 5x
    @test typeof(polynomial(i -> 5, basis)) == polynomialtype(basis, Int)
    @test typeof(polynomial(ones(Int, 1, 1), basis, Int)) <: AbstractPolynomial{Int}
    @test typeof(polynomial(ones(Int, 1, 1), basis, Float64)) <: AbstractPolynomial{Float64}
end

@testset "Enumerate" begin
    monos = [1, x, y^2]
    basis = FixedPolynomialBasis(monos)
    for (i, e) in enumerate(basis) 
        @test e == monos[i]
    end
end

@testset "Coefficients" begin
    @polyvar x y
    p = x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1
    basis = FixedPolynomialBasis([x^4*y^2, x^2*y^4, x^2*y^2, 1])
    coefs = coefficients(p)
    cc = coefficients(p, basis)
    for (i, c) in enumerate(cc)
        @test isapprox(c, coefs[i])
    end

    @test isapprox(p, polynomial(cc, basis ))
    
end
