struct QuotientBasis{T,I,B<:SA.AbstractBasis{T,I},D} <: SA.ExplicitBasis{T,I}
    basis::B
    divisor::D
end

implicit_basis(basis::QuotientBasis) = implicit_basis(basis.basis)
MP.variables(basis::QuotientBasis) = MP.variables(basis.basis)

Base.iterate(basis::QuotientBasis) = iterate(basis.basis)
Base.iterate(basis::QuotientBasis, s) = iterate(basis.basis, s)
Base.length(basis::QuotientBasis) = length(basis.basis)

_object(basis::QuotientBasis) = _object(basis.basis)

function MP.coefficients(p, basis::QuotientBasis)
    return MP.coefficients(rem(p, basis.divisor), basis.basis)
end

function SA.coeffs(coeffs, b::FullBasis{Monomial}, basis::QuotientBasis)
    return MP.coefficients(MP.polynomial(values(coeffs), keys_as_monomials(keys(coeffs), b)), basis)
end

function SA.coeffs(coeffs, sub::SubBasis{Monomial}, basis::QuotientBasis)
    return MP.coefficients(MP.polynomial(coeffs, keys_as_monomials(sub)), basis)
end

function SA.adjoint_coeffs(
    coeffs,
    src::SubBasis{Monomial},
    dest::QuotientBasis{<:Polynomial{Monomial}},
)
    return map(src) do poly
        return sum(MP.terms(rem(MP.polynomial(poly), dest.divisor))) do t
            MP.coefficient(t) * coeffs[SA.key_index(dest.basis, MP.exponents(t))]
        end
    end
end
