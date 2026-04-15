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
    # Converting polynomial to ScaledMonomial and back should roundtrip.
    # This exercises the MStruct because convert_basis uses the algebra
    # multiplication to build the product basis elements.
    p = x^2 + 2x * y + y^2
    a_mono = MB.algebra_element(p)
    scaled_full = MB.FullBasis{MB.ScaledMonomial}(variables(p))
    a_scaled = MB.convert_basis(scaled_full, a_mono)
    # Polynomial roundtrip
    @test polynomial(a_scaled) ≈ p
    # ScaledMonomial coefficients for p = y² + 2xy + x²
    # In ScaledMonomial basis: c(y²)=1, c(√2·xy)=√2, c(x²)=1
    basis_scaled = MB.explicit_basis(a_scaled)
    c_scaled = Vector{Float64}(SA.coeffs(a_scaled, basis_scaled))
    @test c_scaled ≈ [1.0, √2, 1.0]
    # Higher degree: x³y has ScaledMonomial coefficient α/√4 = α/2
    # where scaling(x³y) = √(4!/(3!·1!)) = 2
    q = x^3 * y + x * y^3
    a_q = MB.algebra_element(q)
    a_q_scaled = MB.convert_basis(scaled_full, a_q)
    @test polynomial(a_q_scaled) ≈ q
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
