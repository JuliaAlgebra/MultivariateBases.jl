struct QuotientBasis{T,I,B<:SA.AbstractBasis{T,I},D} <: SA.ExplicitBasis{T,I}
    basis::B
    divisor::D
end

Base.iterate(basis::QuotientBasis) = iterate(basis.basis)
Base.iterate(basis::QuotientBasis, s) = iterate(basis.basis, s)
Base.length(basis::QuotientBasis) = length(basis.basis)

_object(basis::QuotientBasis) = _object(basis.basis)

function MP.coefficients(p, basis::QuotientBasis)
    return MP.coefficients(rem(p, basis.divisor), basis.basis)
end

function SA.coeffs(p, ::FullBasis{Monomial}, basis::Union{QuotientBasis,SubBasis})
    return MP.coefficients(p, basis)
end

function SA.coeffs(coeffs, sub::SubBasis{Monomial}, basis::Union{SubBasis,FullBasis,QuotientBasis})
    return SA.coeffs(MP.polynomial(coeffs, sub.monomials), parent(sub), basis)
end

function MA.operate!(op::SA.UnsafeAddMul{typeof(*)}, a::SA.AlgebraElement{<:Algebra{<:QuotientBasis,Monomial}}, α, p::Polynomial{Monomial})
    _assert_constant(α)
    for t in MP.terms(rem(p.monomial, SA.basis(a).divisor))
        MA.operate!(op, algebra_element(SA.coeffs(a), SA.basis(a).basis), α * MP.coefficient(t), Polynomial{Monomial}(MP.monomial(t)))
        #SA.unsafe_push!(SA.coeffs(a), SA.basis(a).basis[Polynomial{Monomial}(MP.monomial(k))], α * MP.coefficient(t))
    end
    return a
end

function adjoint_coeffs(coeffs, dest::SubBasis{B}, src::SubBasis{B}) where {B<:AbstractMonomial}
    return SA.coeffs(coeffs, dest, src)
end

function adjoint_coeffs(coeffs, dest::QuotientBasis{<:Polynomial{Monomial}}, src::SubBasis{Monomial})
    return map(src.monomials) do mono
        return sum(
            MP.coefficient(t) * coeffs[dest.basis[Polynomial{Monomial}(MP.monomial(t))]]
            for t in MP.terms(rem(mono, dest.divisor))
        )
    end
end
