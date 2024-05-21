using Test
using MultivariateBases
using DynamicPolynomials
@polyvar x y

@testset "Linear" begin
    basis = MB.SubBasis{MB.ScaledMonomial}([x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Float64)
    @test polynomial(i -> i^2, basis) == 4x + y
    @test coefficients(x + 4y, basis) == [4, 1]
end
@testset "Affine" begin
    # It will be sorted and 1 will be moved at the end
    basis = MB.SubBasis{MB.ScaledMonomial}([1, x, y])
    @test polynomial_type(basis, Int) == polynomial_type(x, Float64)
    @test polynomial(i -> i^2, basis) == 9x + 4y + 1
    @test coefficients(9 + x + 4y, basis) == [9, 4, 1]
end
@testset "Quadratic" begin
    basis = MB.SubBasis{MB.ScaledMonomial}([x^2, x * y, y^2])
    @test polynomial_type(basis, Int) == polynomial_type(x, Float64)
    @test polynomial(i -> i^2, basis) == 9x^2 + 4 * √2 * x * y + y^2
    @test coefficients(x^2 + 4x * y + 9y^2, basis) == [9, 4 / √2, 1]
    @test polynomial(basis[1]) == y^2
    @test polynomial(basis[2]) == √2 * x * y
    @test polynomial(basis[3]) == x^2
end
@testset "API degree = $degree" for degree in 0:3
    api_test(MB.ScaledMonomial, degree)
end

@testset "Enumerate" begin
    monos = [1, y, x]
    basis = MB.SubBasis{MB.ScaledMonomial}(monos)
    for (i, e) in enumerate(basis)
        @test polynomial(e) == monos[i]
    end
end

@testset "Coefficients" begin
    coefficient_test(MB.ScaledMonomial, [1, -√3 / √2, 1 / √15, 1 / √15])
end
