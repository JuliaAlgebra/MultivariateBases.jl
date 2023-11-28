using Test
using MultivariateBases
using DynamicPolynomials

@polyvar x y

@testset "Polynomials" begin
    gens = [1 - x^2, x^2 + 2]
    basis = FixedPolynomialBasis(gens)
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(one, basis) == 3
    @test basis[1] == 1 - x^2
    @test basis[2] == x^2 + 2
    @test collect(basis) == gens
    @test generators(basis) == gens
    @test length(basis) == 2
    @test mindegree(basis) == 0
    @test mindegree(basis, x) == 0
    @test maxdegree(basis) == 2
    @test maxdegree(basis, x) == 2
    @test extdegree(basis) == (0, 2)
    @test extdegree(basis, x) == (0, 2)
    @test variables(basis) == [x]
    @test nvariables(basis) == 1
    @test sprint(show, basis) == "FixedPolynomialBasis([1 - x², 2 + x²])"
    @test sprint(print, basis) == "FixedPolynomialBasis([1 - x^2, 2 + x^2])"
    b2 = basis[2:2]
    @test length(b2) == 1
    @test b2[1] == x^2 + 2
    b3 = basis[2:1]
    @test isempty(b3)
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
    @test typeof(polynomial(ones(Int, 1, 1), basis, Int)) <:
          AbstractPolynomial{Int}
    @test typeof(polynomial(ones(Int, 1, 1), basis, Float64)) <:
          AbstractPolynomial{Float64}
end
@testset "Complex" begin
    basis = FixedPolynomialBasis([(1 + 2im) * x])
    @test 5x^2 == @inferred polynomial(ones(Int, 1, 1), basis, Complex{Int})
    @test 5x^2 == @inferred polynomial(ones(Int, 1, 1), basis, Int)
    # TODO not inferred on Julia v1.0
    #@test 5x^2 == @inferred polynomial(ones(Int, 1, 1), basis, Complex{Int})
    #@test 5x^2 == @inferred polynomial(ones(Int, 1, 1), basis, Int)
end
@testset "Empty" begin
    basis = FixedPolynomialBasis(typeof(x + 1)[])
    @test isempty(basis)
    @test isempty(eachindex(basis))
    p = @inferred polynomial(zeros(Int, 0, 0), basis, Int)
    @test iszero(p)
end

@testset "Enumerate" begin
    monos = [1, x, y^2]
    basis = FixedPolynomialBasis(monos)
    for (i, e) in enumerate(basis)
        @test e == monos[i]
    end
end
