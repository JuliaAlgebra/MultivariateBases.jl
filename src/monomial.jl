"""
    struct MonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        monomials::MV
    end

Monomial basis with the monomials of the vector `monomials`.
For instance, `MonomialBasis([1, x, y, x^2, x*y, y^2])` is the monomial basis
for the subspace of quadratic polynomials in the variables `x`, `y`.
"""
struct MonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
    monomials::MV
end
MonomialBasis(monomials::AbstractVector) = MonomialBasis(MP.monovec(monomials))
Base.copy(basis::MonomialBasis) = MonomialBasis(copy(basis.monomials))

Base.length(basis::MonomialBasis) = length(basis.monomials)
empty_basis(::Type{MonomialBasis{MT, MVT}}) where {MT, MVT} = MonomialBasis(MP.emptymonovec(MT))
MP.nvariables(basis::MonomialBasis) = MP.nvariables(basis.monomials)
MP.variables(basis::MonomialBasis) = MP.variables(basis.monomials)
MP.monomialtype(::Type{<:MonomialBasis{MT}}) where MT = MT
MP.polynomialtype(::Union{MonomialBasis{MT}, Type{<:MonomialBasis{MT}}}, T::Type) where MT = MP.polynomialtype(MT, T)
MP.polynomial(f::Function, mb::MonomialBasis) = MP.polynomial(f, mb.monomials)
MP.polynomial(Q::AbstractMatrix, mb::MonomialBasis, T::Type) = MP.polynomial(Q, mb.monomials, T)
function MP.coefficients(p, ::Type{<:MonomialBasis})
    return MP.coefficients(p)
end

# The `i`th index of output is the index of occurence of `x[i]` in `y`,
# or `0` if it does not occur.
function multi_findsorted(x, y)
    I = zeros(Int, length(x))
    j = 1
    for i in eachindex(x)
        while j ≤ length(y) && x[i] < y[j]
            j += 1
        end
        if j ≤ length(y) && x[i] == y[j]
            I[i] = j
        end
    end
    return I
end

function merge_bases(basis1::MonomialBasis, basis2::MonomialBasis)
    monos = MP.mergemonovec([basis1.monomials, basis2.monomials])
    I1 = multi_findsorted(monos, basis1.monomials)
    I2 = multi_findsorted(monos, basis2.monomials)
    return monos, I1, I2
end
