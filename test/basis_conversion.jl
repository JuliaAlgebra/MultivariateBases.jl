module TestBasisConversion

using Test
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import MultivariateBases as MB
import StarAlgebras as SA
import LinearAlgebra

function test_promotion(
    op,
    ::Type{T},
    source::S,
    target::D,
    expected,
) where {T,S,D}
    @test (@inferred SA.promote_coefficient_type(op, T, S, D)) === expected
    SA.promote_coefficient_type(op, T, S, D)
    @test (@allocated SA.promote_coefficient_type(op, T, S, D)) == 0
    return
end

DP.@polyvar x y

@testset "Explicit conversion between basis families" begin
    monomial_basis = MB.FullBasis{MB.Monomial}(x * y)
    for (B, expected) in ((MB.ScaledMonomial, 2x^2), (MB.Chebyshev, 4x^2 - 2))
        full = MB.FullBasis{B}(x * y)
        element = full[[2, 0]]
        t = SA.Term(2, element)
        p = MB.algebra_element(t)
        @test MP.polynomial_type(full, Int) === typeof(p)
        @test MP.polynomial(p) === p
        @test MB.explicit_basis(t) == SA.SubBasis(full, [[2, 0]])
        @test MB.convert_basis(monomial_basis, t) ≈ expected
        @test 2 * MB.convert_basis(monomial_basis, element) ≈ expected
        @test typeof(zero(typeof(p))) === typeof(p)
        @test iszero(zero(typeof(p)))
        for op in (+, -, *)
            @test_throws ArgumentError op(p, MP.polynomial(x))
            @test_throws ArgumentError op(MP.polynomial(x), p)
        end
        sub = MB.SubBasis{B}([x^2, y^2])
        @test MP.variables(sub) == [x, y]
        @test sub.keys == [[0, 2], [2, 0]]
        @test full[[y^2]].keys == [[0, 2]]
    end
end

@testset "Scaled conversion and adjoint" begin
    mono = MB.FullBasis{MB.Monomial}(MP.variables(x * y))
    scaled = MB.FullBasis{MB.ScaledMonomial}(MP.variables(x * y))
    indices = [[0, 0], [0, 1], [1, 0], [1, 1]]
    mono_sub = SA.SubBasis(mono, indices)
    scaled_sub = SA.SubBasis(scaled, indices)
    weights = [1.0, 1.0, 1.0, sqrt(2.0)]
    for (source, target, factors) in (
            (scaled_sub, mono_sub, weights),
            (mono_sub, scaled_sub, inv.(weights)),
        ),
        (T, R) in ((Int, Float64), (BigInt, BigFloat), (Float64, Float64))

        cfs = T[1, 2, 3, 4]
        expected = factors .* cfs
        for op in (SA.coeffs, SA.adjoint_coeffs)
            result = @inferred op(cfs, source, target)
            # Scaling weights are Float64 even for BigInt coefficients.
            @test result ≈ expected rtol = 4eps(Float64)
            @test SA.value_type(result) === R
            test_promotion(op, T, source, target, R)
            empty_result = op(zeros(T, 4), source, target)
            @test iszero(empty_result)
            @test SA.value_type(empty_result) === R
        end
        full_source, full_target = parent(source), parent(target)
        sparse = SA.SparseCoefficients(indices, cfs, SA.comparable(full_source))
        full_result = @inferred SA.coeffs(sparse, full_source, full_target)
        @test [full_result[k] for k in indices] ≈ expected rtol = 4eps(Float64)
        @test SA.value_type(full_result) === R
        @test cfs == T[1, 2, 3, 4]
    end
end

@testset "Reindexing within an alternative basis" begin
    for B in (MB.ScaledMonomial, MB.Chebyshev)
        full = MB.FullBasis{B}(MP.variables(x))
        source = SA.SubBasis(full, [[0], [2]])
        target = SA.SubBasis(full, [[0], [1], [2]])
        cfs = [3, 1]
        @test SA.coeffs(cfs, source, target) == [3, 0, 1]
        @test SA.value_type(SA.coeffs(cfs, source, target)) === Int
        @test SA.adjoint_coeffs([3, 5, 1], source, target) == cfs
        @test SA.value_type(SA.adjoint_coeffs([3, 5, 1], source, target)) ===
              Int
        for op in (SA.coeffs, SA.adjoint_coeffs)
            test_promotion(op, Int, source, target, Int)
            @test op(cfs, source, source) === cfs
            @test op(cfs, source, SA.SubBasis(full, [[0], [2]])) === cfs
        end
        sparse = SA.coeffs(cfs, source, full)
        @test SA.values(sparse) === cfs
        @test SA.keys(sparse) === source.keys
    end
end

@testset "Chebyshev conversion and adjoint" begin
    mono = MB.FullBasis{MB.Monomial}(MP.variables(x))
    cheby = MB.FullBasis{MB.Chebyshev}(MP.variables(x))
    indices = [[0], [1], [2]]
    mono_sub, cheby_sub =
        SA.SubBasis(mono, indices), SA.SubBasis(cheby, indices)
    # T₀ = 1, T₁ = x, T₂ = 2x² - 1.
    A = [1 0 -1; 0 1 0; 0 0 2]
    Ainv = Rational{Int}[1 0 1//2; 0 1 0; 0 0 1//2]
    for (T, R) in
        ((Int, Rational{Int}), (BigInt, Rational{BigInt}), (Float64, Float64))
        cfs = T[1, 2, 3]
        forward = @inferred SA.coeffs(cfs, cheby_sub, mono_sub)
        @test forward == A * cfs
        @test SA.value_type(forward) === R
        reverse = @inferred SA.coeffs(cfs, mono_sub, cheby)
        @test [reverse[k] for k in indices] == Ainv * cfs
        @test SA.value_type(reverse) === R
        for (source, target, transform) in
            ((cheby_sub, mono_sub, A), (mono_sub, cheby_sub, Ainv))
            result = @inferred SA.adjoint_coeffs(cfs, source, target)
            @test result == transform' * cfs
            @test SA.value_type(result) === R
            for op in (SA.coeffs, SA.adjoint_coeffs)
                test_promotion(op, T, source, target, R)
            end
        end
        @test cfs == T[1, 2, 3]
    end
end

end # module
