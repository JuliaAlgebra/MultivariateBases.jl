using Test
using MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Polynomials" begin
    basis = FixedPolynomialBasis([1 - x^2, x^2 + 2])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(one, basis) == 3
end
@testset "Monomial" begin
    basis = FixedPolynomialBasis([x, x^2])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> i^2, basis) == 4x^2 + x
end
@testset "Terms" begin
    basis = FixedPolynomialBasis([1, x, x^2])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> (-1)^i, basis) == -x^2 + x - 1
end
@testset "Linear" begin
    basis = FixedPolynomialBasis([x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(identity, basis) == x + 2y
end
@testset "One variable" begin
    basis = FixedPolynomialBasis([x])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> 5, basis) == 5x
    @test typeof(polynomial(i -> 5, basis)) == polynomial_type(basis, Int)
    @test typeof(polynomial(ones(Int, 1, 1), basis, Int)) <: AbstractPolynomial{Int}
    @test typeof(polynomial(ones(Int, 1, 1), basis, Float64)) <: AbstractPolynomial{Float64}
end
@testset "Complex" begin
    basis = FixedPolynomialBasis([(1 + 2im) * x])
    @test 5x^2 == polynomial(ones(Int, 1, 1), basis, Complex{Int})
    @test 5x^2 == polynomial(ones(Int, 1, 1), basis, Int)
    # TODO not inferred on Julia v1.0
    #@test 5x^2 == @inferred polynomial(ones(Int, 1, 1), basis, Complex{Int})
    #@test 5x^2 == @inferred polynomial(ones(Int, 1, 1), basis, Int)
end
