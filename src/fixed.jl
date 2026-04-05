const FixedBasis = SA.FixedBasis

"""
    struct SimpleBasis{T} <: SA.ExplicitBasis{T,Int}
        elements::Vector{T}
    end

Simple vector-backed basis without star-closure requirement.
Used as inner basis for [`SemisimpleBasis`](@ref).
"""
struct SimpleBasis{T} <: SA.ExplicitBasis{T,Int}
    elements::Vector{T}
end

Base.length(b::SimpleBasis) = length(b.elements)
Base.getindex(b::SimpleBasis, i::Integer) = b.elements[i]
Base.iterate(b::SimpleBasis) = iterate(b.elements)
Base.iterate(b::SimpleBasis, st) = iterate(b.elements, st)

"""
    struct SemisimpleElement{P}
        elements::Vector{P}
    end

Elements of [`SemisimpleBasis`](@ref).
"""
struct SemisimpleElement{P}
    elements::Vector{P}
end
SA.star(p::SemisimpleElement) = SemisimpleElement(SA.star.(p.elements))

function Base.:(==)(a::SemisimpleElement, b::SemisimpleElement)
    return length(a.elements) == length(b.elements) &&
           all(zip(a.elements, b.elements)) do (a, b)
        return a == b
    end
end

function MA.operate!(
    op::SA.UnsafeAddMul,
    res,
    A::SemisimpleElement,
    B::SemisimpleElement,
    α,
)
    for (a, b) in zip(A.elements, B.elements)
        MA.operate!(op, res, a, b, α)
    end
end

"""
    struct SemisimpleBasis{T,B<:SA.ExplicitBasis{T}} <: SA.ExplicitBasis{SemisimpleElement{T},Int}
        bases::Vector{B}
    end

Semisimple basis for use with [SymbolicWedderburn](https://github.com/kalmarek/SymbolicWedderburn.jl/).
Its elements are [`SemisimpleElement`](@ref)s.
"""
struct SemisimpleBasis{T,B<:SA.ExplicitBasis{T}} <: SA.ExplicitBasis{SemisimpleElement{T},Int}
    bases::Vector{B}
end

Base.length(b::SemisimpleBasis) = length(first(b.bases))

function _iterate(b::SemisimpleBasis, elem_state)
    if isnothing(elem_state)
        return
    end
    elem, state = elem_state
    return b[elem], state
end
Base.iterate(b::SemisimpleBasis) = _iterate(b, iterate(Base.OneTo(length(b))))
function Base.iterate(b::SemisimpleBasis, st)
    return _iterate(b, iterate(Base.OneTo(length(b)), st))
end

function Base.getindex(b::SemisimpleBasis, i::Integer)
    return SemisimpleElement(getindex.(b.bases, i))
end

function Base.show(io::IO, b::SemisimpleBasis)
    if length(b.bases) == 1
        print(io, "Simple basis:")
    else
        print(io, "Semisimple basis with $(length(b.bases)) simple sub-bases:")
    end
    for basis in b.bases
        println(io)
        print(io, "  ")
        print(io, basis)
    end
end

# Extract the implicit (full) basis from a SemisimpleBasis
function implicit_basis(
    b::SemisimpleBasis{AE},
) where {AE<:SA.AlgebraElement}
    return implicit_basis(SA.basis(first(first(b.bases))))
end

# parent for SemisimpleBasis - needed for _combine_with_gram
function Base.parent(
    b::SemisimpleBasis{AE},
) where {AE<:SA.AlgebraElement}
    return implicit_basis(b)
end

MP.variables(b::SemisimpleBasis) = MP.variables(implicit_basis(b))
MP.nvariables(b::SemisimpleBasis) = MP.nvariables(implicit_basis(b))
