abstract type AbstractPolynomialVectorBasis{
    PT<:MP.AbstractPolynomialLike,
    PV<:AbstractVector{PT},
} <: AbstractPolynomialBasis end

function MP.monomial_type(
    ::Type{<:AbstractPolynomialVectorBasis{PT}},
) where {PT}
    return MP.monomial_type(PT)
end

function empty_basis(
    B::Type{<:AbstractPolynomialVectorBasis{PT,Vector{PT}}},
) where {PT}
    return B(PT[])
end
function MP.polynomial_type(
    ::AbstractPolynomialVectorBasis{PT},
    T::Type,
) where {PT}
    C = MP.coefficient_type(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomial_type(PT, V)
end
function MP.polynomial(f::Function, basis::AbstractPolynomialVectorBasis{P}) where {P}
    if isempty(generators(basis))
        return zero(P)
    else
        return MP.polynomial(
            mapreduce(
                ip -> f(ip[1]) * ip[2],
                MA.add!!,
                enumerate(basis.polynomials),
            ),
        )
    end
end

function _poly(::MA.Zero, ::Type{P}, ::Type{T}) where {P,T}
    return zero(MP.polynomial_type(P, T))
end

_convert(::Type{P}, p) where {P} = convert(P, p)
_convert(::Type{P}, ::MA.Zero) where {P} = zero(P)

function MP.polynomial(
    Q::AbstractMatrix,
    basis::AbstractPolynomialVectorBasis{P},
    ::Type{T},
) where {P,T}
    n = length(basis)
    @assert size(Q) == (n, n)
    PT = MP.polynomial_type(P, T)
    return _convert(
        PT,
        mapreduce(
            row -> begin
                adjoint(basis.polynomials[row]) * mapreduce(
                    col -> Q[row, col] * basis.polynomials[col],
                    MA.add!!,
                    1:n;
                    init = MA.Zero(),
                   ); end,
            MA.add!!,
            1:n;
            init = MA.Zero(),
        ),
    )::PT
end

"""
    struct FixedPolynomialBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
        polynomials::PV
    end

Polynomial basis with the polynomials of the vector `polynomials`.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct FixedPolynomialBasis{
    PT<:MP.AbstractPolynomialLike,
    PV<:AbstractVector{PT},
} <: AbstractPolynomialVectorBasis{PT,PV}
    polynomials::PV
end
