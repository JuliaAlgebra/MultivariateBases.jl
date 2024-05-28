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

function SA.coeffs(coeffs, q::QuotientBasis, basis::Union{FullBasis, SubBasis})
    # FIXME This is not correct, this is a hack for ConstraintPrimal of SumOfSquares but I should just throw there
    return SA.coeffs(coeffs, q.basis, basis)
end
