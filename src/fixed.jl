"""
    struct FixedBasis{B,M,T,V} <:
        SA.ExplicitBasis{SA.AlgebraElement{Algebra{FullBasis{B,M},B,M},T,V},Int}
        elements::Vector{SA.AlgebraElement{Algebra{FullBasis{B,M},B,M},T,V}}
    end

Fixed basis with polynomials `elements`.
"""
struct FixedBasis{B,M,T,V} <:
       SA.ExplicitBasis{SA.AlgebraElement{Algebra{FullBasis{B,M},B,M},T,V},Int}
    elements::Vector{SA.AlgebraElement{Algebra{FullBasis{B,M},B,M},T,V}}
end

Base.length(b::FixedBasis) = length(b.elements)
Base.getindex(b::FixedBasis, i::Integer) = b.elements[i]

function Base.show(io::IO, b::FixedBasis)
    print(io, "FixedBasis(")
    _show_vector(io, MIME"text/plain"(), b.elements)
    print(io, ")")
    return
end

# TODO refactor with `SA.MappedBasis`
# https://github.com/JuliaAlgebra/StarAlgebras.jl/pull/76

"""
    struct SemisimpleBasis{T,I,B<:SA.ExplicitBasis{T,I}} <: SA.ExplicitBasis{T,I}
        bases::Vector{B}
    end

Semisimple basis for use with [SymbolicWedderburn](https://github.com/kalmarek/SymbolicWedderburn.jl/).
Its elements are of [`SemisimpleElement`](@ref)s.
"""
struct SemisimpleBasis{T,I,B<:SA.ExplicitBasis{T,I}} <: SA.ExplicitBasis{T,I}
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
Base.iterate(b::SemisimpleBasis) = _iterate(b, iterate(keys(first(b.bases))))
function Base.iterate(b::SemisimpleBasis, st)
    return _iterate(b, iterate(keys(first(b.bases)), st))
end

"""
    struct SemisimpleElement{P}
        polynomials::Vector{P}
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
