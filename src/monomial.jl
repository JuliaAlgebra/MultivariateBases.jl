abstract type AbstractMonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis end

Base.length(basis::AbstractMonomialBasis) = length(basis.monomials)
Base.copy(basis::AbstractMonomialBasis) = typeof(basis)(copy(basis.monomials))

MP.nvariables(basis::AbstractMonomialBasis) = MP.nvariables(basis.monomials)
MP.variables(basis::AbstractMonomialBasis) = MP.variables(basis.monomials)
MP.monomialtype(::Type{<:AbstractMonomialBasis{MT}}) where MT = MT

empty_basis(MB::Type{<:AbstractMonomialBasis{MT}}) where {MT} = MB(MP.emptymonovec(MT))
function maxdegree_basis(B::Type{<:AbstractMonomialBasis}, variables, maxdegree::Int)
    return B(MP.monomials(variables, 0:maxdegree))
end
function basis_covering_monomials(B::Type{<:AbstractMonomialBasis}, monos::AbstractVector{<:AbstractMonomial})
    return B(monos)
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

function merge_bases(basis1::MB, basis2::MB) where MB<:AbstractMonomialBasis
    monos = MP.mergemonovec([basis1.monomials, basis2.monomials])
    I1 = multi_findsorted(monos, basis1.monomials)
    I2 = multi_findsorted(monos, basis2.monomials)
    return MB(monos), I1, I2
end

"""
    struct MonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        monomials::MV
    end

Monomial basis with the monomials of the vector `monomials`.
For instance, `MonomialBasis([1, x, y, x^2, x*y, y^2])` is the monomial basis
for the subspace of quadratic polynomials in the variables `x`, `y`.
"""
struct MonomialBasis{MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractMonomialBasis{MT, MV}
    monomials::MV
end
MonomialBasis(monomials::AbstractVector) = MonomialBasis(MP.monovec(monomials))

MP.polynomialtype(::Union{MonomialBasis{MT}, Type{<:MonomialBasis{MT}}}, T::Type) where MT = MP.polynomialtype(MT, T)
MP.polynomial(f::Function, mb::MonomialBasis) = MP.polynomial(f, mb.monomials)

MP.polynomial(Q::AbstractMatrix, mb::MonomialBasis, T::Type) = MP.polynomial(Q, mb.monomials, T)

function MP.coefficients(p, ::Type{<:MonomialBasis})
    return MP.coefficients(p)
end
