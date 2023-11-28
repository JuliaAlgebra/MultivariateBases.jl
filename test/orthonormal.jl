using Test
using MultivariateBases
using DynamicPolynomials

@polyvar x y

@testset "Polynomials" begin
    gens = [1 + x + y + x * y, 1 - x + y - x * y] / 2
    basis = OrthonormalCoefficientsBasis(gens)
    @test iszero(dot(gens[1], gens[2], OrthonormalCoefficientsBasis))
    coefficient_test(basis, [2, -3])
    coefficient_test(basis, [-2im, 1 + 5im])
    coefficient_test(basis, [1im, 2im])
    @test polynomial_type(basis, Int) == polynomial_type(x, Float64)
    @test polynomial(one, basis) == 1 + y
    @test basis[1] == gens[1]
    @test basis[2] == gens[2]
    @test collect(basis) == gens
    @test generators(basis) == gens
    @test length(basis) == 2
    @test mindegree(basis) == 0
    @test mindegree(basis, x) == 0
    @test mindegree(basis, y) == 0
    @test maxdegree(basis) == 2
    @test maxdegree(basis, x) == 1
    @test maxdegree(basis, y) == 1
    @test extdegree(basis) == (0, 2)
    @test extdegree(basis, x) == (0, 1)
    @test extdegree(basis, y) == (0, 1)
    @test variables(basis) == [x, y]
    @test nvariables(basis) == 2
    @test sprint(show, basis) ==
          "OrthonormalCoefficientsBasis([0.5 + 0.5y + 0.5x + 0.5xy, 0.5 + 0.5y - 0.5x - 0.5xy])"
    @test sprint(print, basis) ==
          "OrthonormalCoefficientsBasis([0.5 + 0.5*y + 0.5*x + 0.5*x*y, 0.5 + 0.5*y - 0.5*x - 0.5*x*y])"
    b2 = basis[2:2]
    @test length(b2) == 1
    @test b2[1] == gens[2]
    b3 = basis[2:1]
    @test isempty(b3)
end
@testset "Linear" begin
    basis = OrthonormalCoefficientsBasis([x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(identity, basis) == x + 2y
end
@testset "One variable" begin
    basis = OrthonormalCoefficientsBasis([x])
    @test polynomial_type(basis, Int) == polynomial_type(x, Int)
    @test polynomial(i -> 5, basis) == 5x
    @test typeof(polynomial(i -> 5, basis)) == polynomial_type(basis, Int)
    @test typeof(polynomial(ones(Int, 1, 1), basis, Int)) <:
          AbstractPolynomial{Int}
    @test typeof(polynomial(ones(Int, 1, 1), basis, Float64)) <:
          AbstractPolynomial{Float64}
end
@testset "Complex" begin
    for a in [1, -1, im, -im]
        basis = OrthonormalCoefficientsBasis([a * x])
        @test 5x^2 ==
              @inferred polynomial(5ones(Int, 1, 1), basis, Complex{Int})
        @test 5x^2 == @inferred polynomial(5ones(Int, 1, 1), basis, Int)
        coefficient_test(basis, [2])
        coefficient_test(basis, [-2im])
        coefficient_test(basis, [1 + 5im])
    end
end
@testset "Empty" begin
    basis = OrthonormalCoefficientsBasis(typeof(x + 1)[])
    @test isempty(basis)
    @test isempty(eachindex(basis))
    p = @inferred polynomial(zeros(Int, 0, 0), basis, Int)
    @test iszero(p)
end

@testset "Enumerate" begin
    monos = [1, x, y^2]
    basis = OrthonormalCoefficientsBasis(monos)
    for (i, e) in enumerate(basis)
        @test e == monos[i]
    end
end
