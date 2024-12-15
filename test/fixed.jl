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
    p2 = x_term - im * y_term
    fixed = MB.FixedBasis([p1, p2])
    @test length(fixed) == 2
    @test fixed[1] ≈ p1
    @test fixed[2] ≈ p2
    @test sprint(show, fixed) == "FixedBasis([$p1, $p2])"

    semi = MB.SemisimpleBasis([MB.FixedBasis([p1]), MB.FixedBasis([p2])])
    @test length(semi) == 1
    @test sprint(show, semi) ==
          "Semisimple basis with 2 simple sub-bases:\n  FixedBasis([$p1])\n  FixedBasis([$p2])"
    mult = semi[1]
    @test all(mult.polynomials .≈ [p1, p2])
    smult = SA.star(mult)
    @test all(smult.polynomials .≈ [p2, p1])
end
