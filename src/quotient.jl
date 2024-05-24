struct QuotientBasis{T,I,B<:SA.AbstractBasis{T,I},D} <: SA.ExplicitBasis{T,I}
    basis::B
    divisor::D
end

Base.length(basis::QuotientBasis) = length(basis.basis)

function MP.coefficients(p, basis::QuotientBasis)
    return MP.coefficients(rem(p, basis.divisor), basis.basis)
end
