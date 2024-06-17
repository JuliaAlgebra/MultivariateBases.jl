struct QuotientBasis{T,I,B<:SA.AbstractBasis{T,I},D} <: SA.ExplicitBasis{T,I}
    basis::B
    divisor::D
end

implicit_basis(basis::QuotientBasis) = implicit_basis(basis.basis)

Base.iterate(basis::QuotientBasis) = iterate(basis.basis)
Base.iterate(basis::QuotientBasis, s) = iterate(basis.basis, s)
Base.length(basis::QuotientBasis) = length(basis.basis)

_object(basis::QuotientBasis) = _object(basis.basis)

function MP.coefficients(p, basis::QuotientBasis)
    return MP.coefficients(rem(p, basis.divisor), basis.basis)
end

function SA.coeffs(p, ::FullBasis{Monomial}, basis::QuotientBasis)
    return MP.coefficients(MP.polynomial(p), basis)
end

function SA.coeffs(coeffs, sub::SubBasis{Monomial}, basis::QuotientBasis)
    return SA.coeffs(MP.polynomial(coeffs, sub.monomials), parent(sub), basis)
end

function SA.adjoint_coeffs(
    coeffs,
    src::SubBasis{Monomial},
    dest::QuotientBasis{<:Polynomial{Monomial}},
)
    return map(src.monomials) do mono
        return sum(
            MP.coefficient(t) *
            coeffs[dest.basis[Polynomial{Monomial}(MP.monomial(t))]] for
            t in MP.terms(rem(mono, dest.divisor))
        )
    end
end
