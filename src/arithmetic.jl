# Basis-specific operations that currently expand to the monomial basis.
const _BasisPolynomial{T} =
    SA.AlgebraElement{T,<:SA.StarAlgebra{<:Variables{<:_Basis}}}

function _monomial_expansion(p::_BasisPolynomial)
    return convert_basis(FullBasis{Monomial}(MP.variables(p)), p)
end

function MP.substitute(
    s::MP.AbstractSubstitutionType,
    p::_BasisPolynomial,
    args::MP.Substitutions,
)
    return MP.substitute(s, _monomial_expansion(p), args)
end

function MP.subs(p::_BasisPolynomial, args::MP.AbstractSubstitution...)
    return MP.substitute(MP.Subs(), p, args)
end

function (p::_BasisPolynomial)(args::MP.AbstractSubstitution...)
    return MP.substitute(MP.Eval(), p, args)
end

function (p::_BasisPolynomial)(x::NTuple{N,<:Number}) where {N}
    return (_monomial_expansion(p))(x)
end

function (p::_BasisPolynomial)(x::AbstractVector{<:Number})
    return (_monomial_expansion(p))(x)
end

function (p::_BasisPolynomial)(x::Number...)
    return (_monomial_expansion(p))(x...)
end

function MP.differentiate(p::_BasisPolynomial, x::MP.AbstractVariable)
    return MP.differentiate(_monomial_expansion(p), x)
end
