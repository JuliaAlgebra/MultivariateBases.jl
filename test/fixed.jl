using Test
import StarAlgebras as SA
import MutableArithmetics as MA
using MultivariateBases
const MB = MultivariateBases
using DynamicPolynomials

# Helper: create an implicit AlgebraElement from a coefficient vector and SubBasis,
# copying coefficients to avoid shared-key mutation issues.
function _make_ae(coeffs, basis)
    ae = MB.implicit(MB.algebra_element(collect(coeffs), basis))
    return SA.AlgebraElement(copy(SA.coeffs(ae)), Base.parent(ae))
end

@testset "SimpleBasis" begin
    @polyvar x y
    basis = MB.SubBasis{MB.Monomial}([x^2, x * y, y^2])

    ae1 = _make_ae([1.0, 0.0, 0.0], basis)
    ae2 = _make_ae([0.0, 1.0, 0.0], basis)
    ae3 = _make_ae([0.0, 0.0, 1.0], basis)

    sb = MB.SimpleBasis([ae1, ae2, ae3])

    @testset "length" begin
        @test length(sb) == 3
    end

    @testset "getindex" begin
        @test sb[1] == ae1
        @test sb[2] == ae2
        @test sb[3] == ae3
    end

    @testset "iterate" begin
        elems = collect(sb)
        @test length(elems) == 3
        @test elems[1] == ae1
        @test elems[2] == ae2
        @test elems[3] == ae3
    end

    @testset "single element" begin
        sb1 = MB.SimpleBasis([ae1])
        @test length(sb1) == 1
        @test sb1[1] == ae1
        @test collect(sb1) == [ae1]
    end
end

@testset "SemisimpleElement" begin
    @polyvar x y
    basis = MB.SubBasis{MB.Monomial}([x, y])

    ae1 = _make_ae([1.0, 0.0], basis)
    ae2 = _make_ae([0.0, 1.0], basis)
    ae3 = _make_ae([1.0, 1.0], basis)

    @testset "construction and elements" begin
        se = MB.SemisimpleElement([ae1, ae2])
        @test length(se.elements) == 2
        @test se.elements[1] == ae1
        @test se.elements[2] == ae2
    end

    @testset "equality" begin
        se1 = MB.SemisimpleElement([ae1, ae2])
        se2 = MB.SemisimpleElement([ae1, ae2])
        se3 = MB.SemisimpleElement([ae2, ae1])
        se4 = MB.SemisimpleElement([ae1])
        @test se1 == se2
        @test se1 != se3
        @test se1 != se4
    end

    @testset "star (real)" begin
        se = MB.SemisimpleElement([ae1, ae2])
        s = SA.star(se)
        @test s isa MB.SemisimpleElement
        # For real polynomial bases, star is identity
        @test s == se
    end

    @testset "UnsafeAddMul" begin
        se1 = MB.SemisimpleElement([ae1])
        se2 = MB.SemisimpleElement([ae2])
        alg = MB.algebra(MB.FullBasis{MB.Monomial}(x * y))
        res = zero(Float64, alg)
        MA.operate!(SA.UnsafeAddMul(*), res, se1, se2, 1.0)
        MA.operate!(SA.canonical, res)
        # ae1 ≈ x, ae2 ≈ y, so the product should be x*y
        @test res ≈ x * y
    end

    @testset "UnsafeAddMul with multiple elements" begin
        se1 = MB.SemisimpleElement([ae1, ae3])
        se2 = MB.SemisimpleElement([ae2, ae1])
        alg = MB.algebra(MB.FullBasis{MB.Monomial}(x * y))
        res = zero(Float64, alg)
        MA.operate!(SA.UnsafeAddMul(*), res, se1, se2, 1.0)
        MA.operate!(SA.canonical, res)
        # basis [x,y] is sorted to [y,x], so ae1≈y, ae2≈x, ae3≈y+x
        # sum of products: ae1*ae2 + ae3*ae1 = y*x + (y+x)*y = xy + y^2 + xy = y^2 + 2xy
        @test res ≈ y^2 + 2x * y
    end
end

@testset "SemisimpleBasis" begin
    @polyvar x y
    basis = MB.SubBasis{MB.Monomial}([x^2, x * y, y^2])

    ae1 = _make_ae([1.0, 0.0, 0.0], basis)
    ae2 = _make_ae([0.0, 1.0, 0.0], basis)
    ae3 = _make_ae([0.0, 0.0, 1.0], basis)
    ae4 = _make_ae([1.0, 1.0, 0.0], basis)
    ae5 = _make_ae([0.0, 1.0, 1.0], basis)
    ae6 = _make_ae([1.0, 0.0, 1.0], basis)

    @testset "single simple sub-basis" begin
        sb = MB.SimpleBasis([ae1, ae2, ae3])
        semi = MB.SemisimpleBasis([sb])

        @test length(semi) == 3
        @test length(semi.bases) == 1

        elem = semi[1]
        @test elem isa MB.SemisimpleElement
        @test length(elem.elements) == 1
        @test elem.elements[1] == ae1

        @test semi[2].elements[1] == ae2
        @test semi[3].elements[1] == ae3
    end

    @testset "two simple sub-bases" begin
        sb1 = MB.SimpleBasis([ae1, ae2, ae3])
        sb2 = MB.SimpleBasis([ae4, ae5, ae6])
        semi = MB.SemisimpleBasis([sb1, sb2])

        @test length(semi) == 3
        @test length(semi.bases) == 2

        elem = semi[1]
        @test length(elem.elements) == 2
        @test elem.elements[1] == ae1
        @test elem.elements[2] == ae4

        @test semi[2].elements[1] == ae2
        @test semi[2].elements[2] == ae5

        @test semi[3].elements[1] == ae3
        @test semi[3].elements[2] == ae6
    end

    @testset "iterate" begin
        sb1 = MB.SimpleBasis([ae1, ae2])
        sb2 = MB.SimpleBasis([ae3, ae4])
        semi = MB.SemisimpleBasis([sb1, sb2])

        elems = collect(semi)
        @test length(elems) == 2
        @test elems[1] == semi[1]
        @test elems[2] == semi[2]

        # Test that iteration terminates
        _, st = iterate(semi)
        _, st = iterate(semi, st)
        @test isnothing(iterate(semi, st))
    end

    @testset "show simple" begin
        sb = MB.SimpleBasis([ae1, ae2])
        semi = MB.SemisimpleBasis([sb])
        s = sprint(show, semi)
        @test startswith(s, "Simple basis:")
    end

    @testset "show semisimple" begin
        sb1 = MB.SimpleBasis([ae1, ae2])
        sb2 = MB.SimpleBasis([ae3, ae4])
        semi = MB.SemisimpleBasis([sb1, sb2])
        s = sprint(show, semi)
        @test startswith(s, "Semisimple basis with 2 simple sub-bases:")
    end

    @testset "implicit_basis" begin
        sb = MB.SimpleBasis([ae1, ae2])
        semi = MB.SemisimpleBasis([sb])

        ib = MB.implicit_basis(semi)
        @test ib isa SA.MappedBasis
        @test ib == MB.FullBasis{MB.Monomial}(x * y)
    end

    @testset "parent" begin
        sb = MB.SimpleBasis([ae1, ae2])
        semi = MB.SemisimpleBasis([sb])

        p = parent(semi)
        @test p == MB.implicit_basis(semi)
    end

    @testset "variables" begin
        sb = MB.SimpleBasis([ae1, ae2])
        semi = MB.SemisimpleBasis([sb])

        @test MB.MP.nvariables(semi) == 2
        @test MB.MP.variables(semi) == [x, y]
    end
end

@testset "SemisimpleBasis with complex coefficients" begin
    @polyvar x y
    full = MB.FullBasis{MB.Monomial}(x * y)
    alg = MB.algebra(full)

    # Use 2-variable exponents matching the full basis
    exp_x = MB.MP.exponents(x * y^0)  # exponents for x in x,y variables
    exp_y = MB.MP.exponents(x^0 * y)  # exponents for y in x,y variables

    # Create complex AlgebraElements directly
    ae_cx = SA.AlgebraElement(
        SA.SparseCoefficients(
            [exp_x],
            [1.0 + 1.0im],
            isless,
        ),
        alg,
    )
    ae_cy = SA.AlgebraElement(
        SA.SparseCoefficients(
            [exp_y],
            [1.0 - 1.0im],
            isless,
        ),
        alg,
    )

    sb = MB.SimpleBasis([ae_cx, ae_cy])
    semi = MB.SemisimpleBasis([sb])

    @testset "star of SemisimpleElement with complex" begin
        elem = semi[1]
        s = SA.star(elem)
        @test s isa MB.SemisimpleElement
        # star conjugates: (1+i)*x -> (1-i)*x
        @test s.elements[1] ≈ SA.star(ae_cx)
    end

    @testset "complex element access" begin
        @test length(semi) == 2
        @test semi[1].elements[1] == ae_cx
        @test semi[2].elements[1] == ae_cy
    end

    @testset "complex iterate" begin
        elems = collect(semi)
        @test length(elems) == 2
        @test elems[1].elements[1] == ae_cx
        @test elems[2].elements[1] == ae_cy
    end
end

@testset "FixedBasis alias" begin
    @test MB.FixedBasis === SA.FixedBasis
end
