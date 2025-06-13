struct FullBasis{B<:AbstractMonomialIndexed,M<:MP.AbstractMonomial} <:
       SA.ImplicitBasis{Polynomial{B,M},M} end

function Base.getindex(::FullBasis{B,M}, mono::M) where {B,M}
    return Polynomial{B}(mono)
end

function Base.getindex(::FullBasis{B,M}, p::Polynomial{B,M}) where {B,M}
    return p.monomial
end

# TODO Move it to SA as :
# struct SubBasis{T,I,B<:SA.ImplicitBasis{T,I},V<:AbstractVector{I}} <: SA.ExplicitBasis{T,Int}
#     implicit::B
#     indices::V
# end
struct SubBasis{B<:AbstractMonomialIndexed,M,V<:AbstractVector{M}} <:
       SA.ExplicitBasis{Polynomial{B,M},Int}
    monomials::V
end

# Overload some of the `AbstractVector` interface for convenience
Base.isempty(basis::SubBasis) = isempty(basis.monomials)
Base.eachindex(basis::SubBasis) = eachindex(basis.monomials)
_iterate(::SubBasis, ::Nothing) = nothing
function _iterate(basis::SubBasis{B}, elem_state) where {B}
    return parent(basis)[elem_state[1]], elem_state[2]
end
Base.iterate(basis::SubBasis) = _iterate(basis, iterate(basis.monomials))
Base.iterate(basis::SubBasis, s) = _iterate(basis, iterate(basis.monomials, s))
Base.length(basis::SubBasis) = length(basis.monomials)
Base.firstindex(basis::SubBasis) = firstindex(basis.monomials)
Base.lastindex(basis::SubBasis) = lastindex(basis.monomials)

Base.parent(::SubBasis{B,M}) where {B,M} = FullBasis{B,M}()

function Base.getindex(basis::SubBasis, index::Int)
    return parent(basis)[basis.monomials[index]]
end

function Base.getindex(basis::SubBasis{B,M}, value::Polynomial{B,M}) where {B,M}
    mono = monomial_index(basis, parent(basis)[value])
    if isnothing(mono)
        throw(BoundsError(basis, value))
    end
    return mono
end

const MonomialIndexedBasis{B,M} = Union{SubBasis{B,M},FullBasis{B,M}}
