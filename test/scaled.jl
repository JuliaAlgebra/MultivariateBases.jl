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
    p = x^2 + 4x * y + 9y^2
    @test coefficients(p, basis) == [9, 4 / √2, 1]
    a = MB.algebra_element(p)
    @test SA.coeffs(a, basis) == [9, 4 / √2, 1]
    full = MB.FullBasis{MB.ScaledMonomial}(variables(p))
    @test SA.coeffs(a, full) == SA.SparseCoefficients(
        exponents.(monomial_vector([y^2, x * y, x^2])),
        [9, 4 / √2, 1],
    )
    @test polynomial(basis[1]) == y^2
    @test polynomial(basis[2]) == √2 * x * y
    @test polynomial(basis[3]) == x^2
end
@testset "coeffs(Polynomial{ScaledMonomial}, FullBasis{Monomial})" begin
    mono_full = MB.FullBasis{MB.Monomial}([x, y])
    # ScaledMonomial(xy) = √2·xy, so coefficient in Monomial basis is √2
    p_xy = collect(MB.SubBasis{MB.ScaledMonomial}([x * y]))[1]
    c_xy = SA.coeffs(p_xy, mono_full)
    @test only(SA.nonzero_pairs(c_xy))[2] ≈ √2
    # ScaledMonomial(x²) = x², so coefficient in Monomial basis is 1
    p_x2 = collect(MB.SubBasis{MB.ScaledMonomial}([x^2]))[1]
    c_x2 = SA.coeffs(p_x2, mono_full)
    @test only(SA.nonzero_pairs(c_x2))[2] ≈ 1.0
    # ScaledMonomial(x²y) = √3·x²y, coefficient is √3
    p_x2y = collect(MB.SubBasis{MB.ScaledMonomial}([x^2 * y]))[1]
    c_x2y = SA.coeffs(p_x2y, mono_full)
    @test only(SA.nonzero_pairs(c_x2y))[2] ≈ √3
end
@testset "Algebra multiplication (MStruct)" begin
    # The MStruct for ScaledMonomial computes the structure coefficient
    # c = scaling(a)·scaling(b) / scaling(a+b) for the product ŝ_a · ŝ_b = c · ŝ_{a+b}.
    # Test this directly by multiplying algebra elements and inspecting coefficients.
    scaled_full = MB.FullBasis{MB.ScaledMonomial}(variables(x * y))
    alg = MB.algebra(scaled_full)
    e_y = SA.AlgebraElement(SA.SparseCoefficients([[0, 1]], [1.0]), alg)
    e_x = SA.AlgebraElement(SA.SparseCoefficients([[1, 0]], [1.0]), alg)
    # ŷ · x̂: scaling(y)=1, scaling(x)=1, scaling(xy)=√2 → coeff = 1/√2
    prod_yx = e_y * e_x
    pairs_yx = collect(SA.nonzero_pairs(SA.coeffs(prod_yx)))
    @test length(pairs_yx) == 1
    @test pairs_yx[1][1] == [1, 1]
    @test pairs_yx[1][2] ≈ 1 / √2
    # ŷ · ŷ: scaling(y)=1, scaling(y²)=1 → coeff = 1
    prod_yy = e_y * e_y
    pairs_yy = collect(SA.nonzero_pairs(SA.coeffs(prod_yy)))
    @test length(pairs_yy) == 1
    @test pairs_yy[1][1] == [0, 2]
    @test pairs_yy[1][2] ≈ 1.0
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
