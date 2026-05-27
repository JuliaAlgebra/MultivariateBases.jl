using Test
import StarAlgebras as SA
import MutableArithmetics as MA
using MultivariateBases
const MB = MultivariateBases

@testset "implicit_basis on generic SA.SubBasis" begin
    # `MB.implicit_basis` was historically only defined on `MB.SubBasis`
    # (= `SA.SubBasis{Polynomial}`). It is now defined on any `SA.SubBasis`,
    # so that a `SubBasis` whose parent is a non-polynomial `SA.ImplicitBasis`
    # (such as the `SA.DiracBasis` wrapping a custom algebra of monoid elements
    # used by the CHSH tutorial) routes through the same fallback.
    parent_basis = SA.DiracBasis(["a", "b", "c"])
    sub = SA.SubBasis(parent_basis, ["a", "b"])

    @test sub isa SA.SubBasis
    @test !(sub isa MB.SubBasis)  # not a polynomial-typed SubBasis
    @test MB.implicit_basis(sub) === parent_basis
    @test MA.promote_operation(MB.implicit_basis, typeof(sub)) ==
          typeof(parent_basis)
end
