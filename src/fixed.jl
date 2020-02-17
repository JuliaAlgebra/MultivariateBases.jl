abstract type AbstractPolynomialVectorBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis end

Base.length(basis::AbstractPolynomialVectorBasis) = length(basis.polynomials)
Base.copy(basis::AbstractPolynomialVectorBasis) = typeof(basis)(copy(basis.polynomials))

MP.nvariables(basis::AbstractPolynomialVectorBasis) = MP.nvariables(basis.polynomials)
MP.variables(basis::AbstractPolynomialVectorBasis) = MP.variables(basis.polynomials)
MP.monomialtype(::Type{<:AbstractPolynomialVectorBasis{PT}}) where PT = MP.monomialtype(PT)

empty_basis(B::Type{<:AbstractPolynomialVectorBasis{PT, Vector{PT}}}) where PT = B(PT[])
function MP.polynomialtype(basis::AbstractPolynomialVectorBasis{PT}, T::Type) where PT
    C = MP.coefficienttype(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomialtype(PT, V)
end
function MP.polynomial(f::Function, basis::AbstractPolynomialVectorBasis)
    return MP.polynomial(mapreduce(
        ip -> f(ip[1]) * ip[2], MA.add!, enumerate(basis.polynomials)))
end

function MP.polynomial(Q::AbstractMatrix, basis::AbstractPolynomialVectorBasis,
                       T::Type)
    n = length(basis)
    @assert size(Q) == (n, n)
    return MP.polynomial(mapreduce(row -> basis.polynomials[row] *
        mapreduce(col -> Q[row, col] * basis.polynomials[col], MA.add!, 1:n),
        MA.add!, 1:n), T)
end

"""
    struct FixedPolynomialBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
        polynomials::PV
    end

Polynomial basis with the polynomials of the vector `polynomials`.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct FixedPolynomialBasis{PT<:MP.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialVectorBasis{PT, PV}
    polynomials::PV
end
