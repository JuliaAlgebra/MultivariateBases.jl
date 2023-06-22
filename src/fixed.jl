abstract type AbstractPolynomialVectorBasis{
    PT<:MP.AbstractPolynomialLike,
    PV<:AbstractVector{PT},
} <: AbstractPolynomialBasis end

Base.length(basis::AbstractPolynomialVectorBasis) = length(basis.polynomials)
function Base.copy(basis::AbstractPolynomialVectorBasis)
    return typeof(basis)(copy(basis.polynomials))
end

Base.firstindex(basis::AbstractPolynomialVectorBasis) = 1
Base.lastindex(basis::AbstractPolynomialVectorBasis) = length(basis)
function Base.getindex(basis::AbstractPolynomialVectorBasis, i::Int)
    return basis.polynomials[i]
end

function MP.nvariables(basis::AbstractPolynomialVectorBasis)
    return MP.nvariables(basis.polynomials)
end
function MP.variables(basis::AbstractPolynomialVectorBasis)
    return MP.variables(basis.polynomials)
end
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
function MP.polynomial(f::Function, basis::AbstractPolynomialVectorBasis)
    return MP.polynomial(
        mapreduce(
            ip -> f(ip[1]) * ip[2],
            MA.add!!,
            enumerate(basis.polynomials),
        ),
    )
end

function MP.polynomial(
    Q::AbstractMatrix,
    basis::AbstractPolynomialVectorBasis,
    T::Type,
)
    n = length(basis)
    @assert size(Q) == (n, n)
    return MP.polynomial(
        mapreduce(
            row ->
                adjoint(basis.polynomials[row]) * mapreduce(
                    col -> Q[row, col] * basis.polynomials[col],
                    MA.add!!,
                    1:n,
                ),
            MA.add!!,
            1:n,
        ),
        T,
    )
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
