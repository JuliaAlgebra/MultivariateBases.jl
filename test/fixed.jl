using Test
import StarAlgebras as SA
using MultivariateBases
const MB = MultivariateBases
using DynamicPolynomials

@testset "FixedBasis" begin
    @polyvar x y
    x_term = MB.term_element(1, MB.Polynomial{MB.Monomial}(x))
    y_term = MB.term_element(1, MB.Polynomial{MB.Monomial}(y))
    p1 = x_term + im * y_term
    q1 = x_term - im * y_term
    p2 = im * x_term - 2y_term
    q2 = -im * x_term - 2y_term
    fixed = MB.FixedBasis([p1, p2])
    @test length(fixed) == 2
    @test fixed[1] ≈ p1
    @test fixed[2] ≈ p2
    @test sprint(show, fixed) == "FixedBasis([$p1, $p2])"

    semi = MB.SemisimpleBasis([MB.FixedBasis([p1, p2]), MB.FixedBasis([q1, q2])])
    @test length(semi) == 2
    @test sprint(show, semi) ==
          "Semisimple basis with 2 simple sub-bases:\n  FixedBasis([$p1, $p2])\n  FixedBasis([$q1, $q2])"
    mult = semi[1]
    @test all(mult.elements .≈ [p1, q1])
    smult = SA.star(mult)
    @test all(smult.elements .≈ [q1, p1])
    res = zero(Complex{Int}, MB.algebra(MB.FullBasis{Monomial,typeof(x*y)}()));
    MA.operate!(SA.UnsafeAddMul(*), res, semi[1], semi[2], true)
    MA.operate!(SA.canonical, res)
    @test res == p1 * p2 + q1 * q2
    MA.operate!(SA.UnsafeAddMul(*), res, semi[1], semi[2], false)
    @test res == p1 * p2 + q1 * q2
    MA.operate!(SA.UnsafeAddMul(*), res, semi[2], semi[1], -1)
    MA.operate!(SA.canonical, res)
    @test iszero(res)
    MA.operate!(SA.UnsafeAddMul(*), res, semi[2], semi[2], -1)
    MA.operate!(SA.canonical, res)
    @test res == -p2^2 - q2^2
end
