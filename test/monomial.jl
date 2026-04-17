using Test
using MultivariateBases
const MB = MultivariateBases
using MultivariatePolynomials

function test_monomial(x, y)
    @testset "StarAlgebras" begin
        vars = MB.Variables{MB.Monomial}(variables(x))
        @test vars(exponents(x^2)) == vars(exponents(x^2))
        @test vars(exponents(x^3)) != vars(exponents(x^2))
        o = vars(exponents(constant_monomial(x^2)))
        @test isone(o)
    end

    @testset "Linear" begin
        basis = SubBasis{MB.Monomial}([x, y])
        @test basis == SubBasis{MB.Monomial}([y, x])
        @test basis != SubBasis{MB.Monomial}([y, y^2, x])
        @test polynomial_type(basis, Int) == polynomial_type(x * y, Int)
        @test polynomial(i -> i^2, basis) == 4x + y
        @test coefficients(x + 4y, basis) == [4, 1]
        @test polynomial(basis[1]) == y
        @test polynomial(basis[2]) == x
        @test MB.keys_as_monomials(basis) == [y, x]
        @test polynomial.(collect(basis)) == [y, x]
        @test variables(basis) == variables(x * y)
        @test nvariables(basis) == 2
        @test sprint(show, basis) == "SubBasis{Monomial}([y, x])"
        @test sprint(show, MIME"text/print"(), basis) ==
              "SubBasis{Monomial}([y, x])"
        @test sprint(show, MIME"text/plain"(), basis) ==
              "SubBasis{Monomial}([y, x])"
        @test sprint(print, basis) == "SubBasis{Monomial}([y, x])"
    end

    @testset "FullBasis type stability" begin
        fb = @inferred MB.FullBasis{MB.Monomial}(x * y)
        @test isconcretetype(typeof(fb))
        a = @inferred MB.algebra_element(x + y)
        @test isconcretetype(typeof(a))
    end

    @testset "FullBasis printing" begin
        fb = MB.FullBasis{MB.Monomial}(x * y)
        @test sprint(show, fb) == "FullBasis{Monomial}([x, y])"
        @test sprint(show, MIME"text/print"(), fb) ==
              "FullBasis{Monomial}([x, y])"
        @test sprint(show, MIME"text/plain"(), fb) ==
              "FullBasis{Monomial}([x, y])"
        @test sprint(print, fb) == "FullBasis{Monomial}([x, y])"
        # Single variable
        fb1 = MB.FullBasis{MB.Monomial}(x)
        @test sprint(show, fb1) == "FullBasis{Monomial}([x])"
        @test sprint(print, fb1) == "FullBasis{Monomial}([x])"
    end

    @testset "Affine" begin
        # It will be sorted and 1 will be moved at the end
        basis = SubBasis{MB.Monomial}([1, x, y])
        @test polynomial_type(basis, Int) == polynomial_type(x * y, Int)
        poly = 9x + 4y + 1
        @test polynomial(i -> i^2, basis) == poly
        @test coefficients(9 + x + 4y, basis) == [9, 4, 1]
        @test MB.explicit_basis(poly) == basis
    end

    @testset "Quadratic" begin
        basis = SubBasis{MB.Monomial}([x^2, x * y, y^2])
        @test polynomial_type(basis, Int) == polynomial_type(x * y, Int)
        @test polynomial(i -> i^2, basis) == 9x^2 + 4x * y + y^2
        @test coefficients(x^2 + 4x * y + 9y^2, basis) == [9, 4, 1]
        @test sprint(show, basis) == "SubBasis{Monomial}([y², xy, x²])"
        @test sprint(print, basis) == "SubBasis{Monomial}([y^2, x*y, x^2])"
    end

    @testset "API degree = $degree" for degree in 0:3
        api_test(MB.Monomial, degree)
    end

    @testset "Empty" begin
        basis = MB.SubBasis{MB.Monomial}(typeof(x^2)[])
        @test isempty(basis)
        @test isempty(eachindex(basis))
        p = @inferred polynomial(zeros(Int, 0, 0), basis, Int)
        @test iszero(p)
    end

    @testset "Enumerate" begin
        monos = [1, y, x]
        basis = MB.SubBasis{MB.Monomial}(monos)
        for (i, e) in enumerate(basis)
            @test polynomial(e) == monos[i]
        end
    end

    @testset "Coefficients" begin
        coefficient_test(MB.Monomial, [1, -3, 1, 1])
    end

    # All SparseCoefficients in MultivariateBases must use the algebra's
    # Graded{LexOrder}, not Julia's default `isless` on exponent vectors.
    # `isless` gives lex order [0,2] < [1,0], but Graded{LexOrder} gives
    # [1,0] < [0,2] (lower total degree first). These disagree when mixing
    # monomials of different total degrees (e.g. x degree 1 vs y² degree 2).
    # The helper below checks that after `canonical`, keys stay in the
    # algebra's ordering. With the old `isless` default, canonical would
    # reorder [1,0],[0,2] to [0,2],[1,0], silently swapping coefficients.
    function test_graded_order(ae)
        MA.operate!(SA.canonical, SA.coeffs(ae))
        keys = collect(SA.keys(SA.coeffs(ae)))
        vals = collect(SA.values(SA.coeffs(ae)))
        return keys, vals
    end

    @testset "sparse_coefficients ordering" begin
        p = x + y^2  # mixes degree 1 (x) and degree 2 (y²)
        full = MB.FullBasis{MB.Monomial}([x, y])
        ae = MB.algebra_element(MB.sparse_coefficients(p), full)
        keys, vals = test_graded_order(ae)
        @test keys == [[1, 0], [0, 2]]
        @test vals == [1, 1]
        eb = MB.explicit_basis(ae)
        @test collect(eb.keys) == [[1, 0], [0, 2]]
    end

    @testset "term_element ordering" begin
        # term_element is used when multiplying a scalar by a Polynomial
        full = MB.FullBasis{MB.Monomial}([x, y])
        p_x = full[[1, 0]]    # the Polynomial for x
        ae = 2 * p_x           # calls term_element
        keys, _ = test_graded_order(ae)
        @test keys == [[1, 0]]
        # Now test with y² to make sure the ordering type is right
        p_y2 = full[[0, 2]]   # the Polynomial for y²
        # Combine: this triggers canonical on mixed-degree keys
        combined = ae + (3 * p_y2)
        keys, vals = test_graded_order(combined)
        @test keys == [[1, 0], [0, 2]]
        @test vals == [2, 3]
    end

    @testset "constant_algebra_element ordering" begin
        full = MB.FullBasis{MB.Monomial}([x, y])
        ce = MB.constant_algebra_element(full, Float64)
        keys, vals = test_graded_order(ce)
        @test keys == [[0, 0]]
        @test vals == [1.0]
        # Combine with a degree-2 element to trigger ordering
        p_y2 = full[[0, 2]]
        combined = ce + (2.0 * p_y2)
        keys, vals = test_graded_order(combined)
        @test keys == [[0, 0], [0, 2]]
        @test vals == [1.0, 2.0]
    end

    @testset "MStruct multiplication ordering" begin
        full = MB.FullBasis{MB.Monomial}([x, y])
        # (xy + 1)(x + y²) = x²y + xy³ + x + y²
        # Degrees 3,4,1,2 → Graded: [1,0],[0,2],[2,1],[1,3]
        a = MB.algebra_element(MB.sparse_coefficients(x * y + 1 + 0 * x), full)
        b = MB.algebra_element(MB.sparse_coefficients(x + y^2), full)
        c = a * b
        keys, vals = test_graded_order(c)
        @test keys == [[1, 0], [0, 2], [2, 1], [1, 3]]
        @test vals == [1, 1, 1, 1]
    end

    @testset "coeffs SubBasis→FullBasis ordering" begin
        full = MB.FullBasis{MB.Monomial}([x, y])
        sub = MB.SubBasis{MB.Monomial}([x, y^2])  # mixed degrees
        cfs = SA.coeffs([3, 7], sub, full)
        MA.operate!(SA.canonical, cfs)
        keys = collect(SA.keys(cfs))
        vals = collect(SA.values(cfs))
        @test keys == [[1, 0], [0, 2]]
        @test vals == [3, 7]
    end

    @testset "promote_bases_with_maps" begin
        a = MB.algebra_element(x - x^2)
        b = MB.FullBasis{MB.Monomial}(x * y)
        a2, b2 = SA.promote_bases(a, b)
        @test collect.(a2.coeffs.basis_elements) == [[1, 0], [2, 0]]
        @test collect.(a2.coeffs.basis_elements) == [[1, 0], [2, 0]]
        @test SA.basis(a) != b
        @test SA.basis(a2) == b
        b2 === b
    end

    @testset "promote_bases_with_maps MStruct" begin
        # Test MStruct-MStruct promotion with different variables
        a = MB.algebra_element(x - x^2)
        b = MB.algebra_element(y + y^2)
        alg_a = SA.parent(a)
        alg_b = SA.parent(b)
        mstruct_a = alg_a.mstructure
        mstruct_b = alg_b.mstructure
        @test mstruct_a isa MB.MStruct
        @test mstruct_b isa MB.MStruct
        (new_a, map_a), (new_b, map_b) =
            SA.promote_bases_with_maps(mstruct_a, mstruct_b)
        # After promotion, both MStructs should have a basis
        # with variables [x, y]
        @test new_a isa MB.MStruct
        @test new_b isa MB.MStruct
        @test variables(SA.basis(new_a)) == variables(x * y)
        @test variables(SA.basis(new_b)) == variables(x * y)
        # Test that maps are callable and remap indices correctly
        @test map_a isa Function
        @test map_b isa Function
        # Test MStruct-MStruct promotion with same variables (identity case)
        c = MB.algebra_element(x + y)
        alg_c = SA.parent(c)
        mstruct_c = alg_c.mstructure
        (new_c1, map_c1), (new_c2, map_c2) =
            SA.promote_bases_with_maps(mstruct_c, mstruct_c)
        @test variables(SA.basis(new_c1)) == variables(x * y)
        @test variables(SA.basis(new_c2)) == variables(x * y)
        # Identity promotion should preserve the basis
        @test SA.basis(new_c1) == SA.basis(mstruct_c)
        @test SA.basis(new_c2) == SA.basis(mstruct_c)
    end

    @testset "promote_bases MStruct across variables" begin
        # Test that promote_bases works end-to-end on algebra elements
        # with disjoint variables (exercises MStruct-MStruct promotion)
        a = MB.algebra_element(x - x^2)
        b = MB.algebra_element(y + y^2)
        a2, b2 = SA.promote_bases(a, b)
        @test variables(a2) == variables(x * y)
        @test variables(b2) == variables(x * y)
        @test polynomial(a2) ≈ x - x^2
        @test polynomial(b2) ≈ y + y^2
        # After promotion, parents should match
        @test SA.parent(a2) == SA.parent(b2)
        # Test with overlapping but non-identical variables
        c = MB.algebra_element(x + y)
        d = MB.algebra_element(x^2)
        c2, d2 = SA.promote_bases(c, d)
        @test variables(c2) == variables(x * y)
        @test variables(d2) == variables(x * y)
        @test polynomial(c2) ≈ x + y
        @test polynomial(d2) ≈ x^2
        @test SA.parent(c2) == SA.parent(d2)
        # Promoted algebra elements can do arithmetic in the same algebra
        s = a2 + b2
        @test polynomial(s) ≈ x - x^2 + y + y^2
    end

    @testset "hash" begin
        monosx = [1, x]
        basisx = MB.SubBasis{MB.Monomial}(monosx)
        scaled_basisx = MB.SubBasis{MB.ScaledMonomial}(monosx)
        monosy = [1, y]
        basisy = MB.SubBasis{MB.Monomial}(monosy)
        @test hash(basisx[2]) != hash(scaled_basisx[2])
        @test hash(basisx[2]) != hash(basisy[2])
    end
end

import DynamicPolynomials
test_monomial((DynamicPolynomials.@polyvar x y)...)

import TypedPolynomials
test_monomial((TypedPolynomials.@polyvar x y)...)

import MutableArithmetics as MA

@testset "Noncommutative" begin
    DynamicPolynomials.@ncpolyvar x y

    @testset "NC FullBasis" begin
        fb = MB.FullBasis{MB.Monomial}(x * y)
        @test variables(fb) == [x, y]
        # NC FullBasis uses DiracMStructure, not MStruct
        alg = MB.algebra(fb)
        @test alg.mstructure isa SA.DiracMStructure
    end

    @testset "NC SubBasis" begin
        sb = MB.SubBasis{MB.Monomial}([x * y, y * x, x^2])
        @test length(sb) == 3
        monos = MB.keys_as_monomials(sb)
        # xy and yx are distinct NC monomials
        @test x * y in monos
        @test y * x in monos
        @test x^2 in monos
    end

    @testset "NC algebra_element round-trip" begin
        p = x * y + y * x + 2x^2
        a = MB.algebra_element(p)
        # Round-trip: polynomial → algebra_element → polynomial
        @test MultivariatePolynomials.polynomial(a) == p
        # NC terms are kept separate: xy, yx, x^2
        @test length(SA.supp(a)) == 3
    end

    @testset "NC monomial distinctness" begin
        # In NC, xy ≠ yx
        @test x * y != y * x
        a_xy = MB.algebra_element(x * y)
        a_yx = MB.algebra_element(y * x)
        @test !(a_xy ≈ a_yx)
    end

    @testset "NC degree-2 polynomials" begin
        # Various degree-2 NC polynomials preserve structure
        for p in
            [x^2 + y^2, x * y + y * x, x * y - y * x, x^2 + x * y + y * x + y^2]
            a = MB.algebra_element(p)
            @test MultivariatePolynomials.polynomial(a) == p
        end
    end

    @testset "NC explicit_basis" begin
        p = x * y + y * x + x^2 + y^2
        a = MB.algebra_element(p)
        eb = MB.explicit_basis(a)
        @test length(eb) == 4
        monos = MB.keys_as_monomials(eb)
        @test x * y in monos
        @test y * x in monos
        @test x^2 in monos
        @test y^2 in monos
    end

    @testset "NC Polynomial type" begin
        p_xy = MB.Polynomial{MB.Monomial}(x * y)
        p_yx = MB.Polynomial{MB.Monomial}(y * x)
        # These should be distinct NC Polynomial wrappers
        @test p_xy != p_yx
        @test monomial(p_xy) == x * y
        @test monomial(p_yx) == y * x
    end

    @testset "NC algebra_element arithmetic" begin
        # Scalar multiplication preserves NC structure
        a = MB.algebra_element(x * y)
        c = 3 * a
        @test MultivariatePolynomials.polynomial(c) == 3 * x * y
        # Addition/subtraction of NC polynomials (same variables, same algebra)
        p = x * y + y * x
        a_sum = MB.algebra_element(p)
        @test MultivariatePolynomials.polynomial(2 * a_sum) == 2p
        @test MultivariatePolynomials.polynomial(a_sum - a_sum) == 0 * p
    end
end
