# Inspired by `Hypatia.jl/src/PolyUtils/realinterp.jl`

export BoxSampling

abstract type AbstractNodes{T,V} end

mutable struct BoxSampling{T,V} <: AbstractNodes{T,V}
    lower::V
    upper::V
    sample_factor::Int
    orthogonalize::Bool
    function BoxSampling(lower, upper; sample_factor = 0, orthogonalize = false)
        @assert length(lower) == length(upper)
        l = float.(lower)
        u = float.(upper)
        V = promote_type(typeof(l), typeof(u))
        return new{eltype(V),V}(l, u, sample_factor, orthogonalize)
    end
end

function sample(s::BoxSampling, n::Integer) where {T}
    samples = rand(T, length(s.lower), n) .- inv(T(2))
    shift = (dom.u .+ dom.l) .* inv(T(2))
    for i in 1:n
        for j in eachindex(s.lower)
            samples[j, i] = samples[j, i] * (dom.u[j] - dom.l[j]) + shift[j]
        end
    end
    return samples
end

struct LagrangePolynomial{T,V}
    point::V
end

struct ImplicitLagrangeBasis{T,V,N<:AbstractNodes{T,V}} <:
       SA.ImplicitBasis{LagrangePolynomial{T,V},V}
    sampling::AbstractNodes{T,V}
    function ImplicitLagrangeBasis(nodes::AbstractNodes{T,V}) where {T,V}
        return new{T,V,typeof(nodes)}(nodes)
    end
end
struct LagrangeBasis{T,P,V<:AbstractVector{P}} <:
       SA.ExplicitBasis{LagrangePolynomial{T,V},Int}
    points::V
end

Base.length(basis::LagrangeBasis) = length(basis.points)

function transformation_to(basis::SubBasis, lag::LagrangeBasis{T}) where {T}
    V = Matrix{T}(undef, length(lag), length(basis))
    for 
end

function sample(s::AbstractNodes, basis::SubBasis)
    samples = sample(s, length(basis) * s.sample_factor)
    return eachcol(samples)
end

function explicit_basis_covering(
    implicit::ImplicitLagrangeBasis,
    basis::SubBasis,
)
    return LagrangeBasis(sample(implicit.sampling, length(basis)))
end
