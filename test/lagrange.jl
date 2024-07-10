module TestLagrange

using Random, Test
import StarAlgebras as SA
using DynamicPolynomials
import MultivariateBases as MB

function _test(B::Type)
    @polyvar x[1:2]
    Random.seed!(0)
    implicit =
        MB.ImplicitLagrangeBasis(x, MB.BoxSampling([-1, -1], UInt32[1, 1]))
    point = zeros(2)
    poly = implicit[x=>point]
    @test poly isa MB.LagrangePolynomial
    @test poly.variables == x
    @test poly.point === point
    err = ErrorException(
        "Variables `$([x[1]])` do not match Lagrange basis variables `$x`",
    )
    @test_throws err implicit[[x[1]]=>[0.0]]
    monos = monomials(x, 0:2)
    coeffs = collect(eachindex(monos))
    sub = MB.SubBasis{B}(monos)
    lag = MB.explicit_basis_covering(implicit, sub)
    @test typeof(lag) == MB.explicit_basis_type(typeof(implicit))
    a = MB.algebra_element(
        SA.SparseCoefficients(collect(monos), coeffs),
        MB.FullBasis{B,typeof(prod(x))}(),
    )
    @test SA.coeffs(a, lag) ≈ SA.coeffs(coeffs, sub, lag)
    @polyvar z
    bad = MB.SubBasis{B}([prod(x) * z])
    err = ErrorException(
        "Cannot evaluate `$bad` as its variable `$z` is not part of the variables `$x` of the `LagrangeBasis`",
    )
    @test_throws err SA.coeffs([1], bad, lag)
end

test_mono() = _test(MB.Monomial)

test_cheby() = _test(MB.Chebyshev)

function test_num_samples()
    @test MB.num_samples(1, 1) == 1
    @test MB.num_samples(0, 1) == 10
    @test MB.num_samples(0, 15_000) == 5 * 15_000
    @test MB.num_samples(0, 22_000) == 2 * 22_000
    @test MB.num_samples(0, 23_000) == 23_000
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestLagrange.runtests()
