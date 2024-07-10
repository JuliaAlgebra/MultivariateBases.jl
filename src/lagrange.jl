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

function sample(s::BoxSampling{T}, n::Integer) where {T}
    samples = rand(T, length(s.lower), n) .- inv(T(2))
    shift = (s.upper .+ s.lower) .* inv(T(2))
    for i in 1:n
        for j in eachindex(s.lower)
            samples[j, i] = samples[j, i] * (s.upper[j] - s.lower[j]) + shift[j]
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
    function LagrangeBasis(points::AbstractVector)
        P = eltype(points)
        return new{eltype(P), P, typeof(points)}(points)
    end
end

Base.length(basis::LagrangeBasis) = length(basis.points)
function Base.getindex(basis::LagrangeBasis, I::AbstractVector{<:Integer})
    return LagrangeBasis(basis.points[I])
end

function eval_basis!(univariate_buffer, result, basis::SubBasis{B}, values) where {B}
    for i in eachindex(values)
        univariate_eval!(B, view(univariate_buffer, :, i), values[i])
    end
    for i in eachindex(basis)
        result[i] = one(eltype(result))
        exp = MP.exponents(basis.monomials[i])
        @assert length(exp) == length(values)
        for j in eachindex(values)
            result[i] = MA.operate!!(*, result[i], univariate_buffer[exp[j] + 1, j])
        end
    end
    return result
end

function transformation_to(basis::SubBasis, lag::LagrangeBasis{T}) where {T}
    # To avoid allocating this too often, we allocate it once here
    # and reuse it for each sample
    univariate_buffer = Matrix{T}(undef, length(basis), MP.nvariables(basis))
    V = Matrix{T}(undef, length(lag), length(basis))
    for i in eachindex(lag)
        eval_basis!(univariate_buffer, view(V, i, :), basis, lag.points[i])
    end
    return V
end

# Heuristic taken from Hypatia
function num_samples(sample_factor, dim)
    if iszero(sample_factor)
        if dim <= 12_000
            sample_factor = 10
        elseif dim <= 15_000
            sample_factor = 5
        elseif dim <= 22_000
            sample_factor = 2
        else
            sample_factor = 1
        end
    end
    return sample_factor * dim
end

function sample(s::AbstractNodes, basis::SubBasis)
    full = LagrangeBasis(eachcol(sample(s, num_samples(s.sample_factor, length(basis)))))
    V = transformation_to(basis, full)
    display(V)
    F = LinearAlgebra.qr!(Matrix(V'), LinearAlgebra.ColumnNorm())
    display(F)
    kept_indices = F.p[1:length(basis)]
    return full[kept_indices]
end

function explicit_basis_covering(
    implicit::ImplicitLagrangeBasis,
    basis::SubBasis,
)
    return sample(implicit.sampling, basis)
end

function SA.coeffs(coeffs, source::SubBasis, target::LagrangeBasis)
    return transformation_to(source, target) * coeffs
end

function SA.coeffs(coeffs, implicit::FullBasis, target::LagrangeBasis)
    a = algebra_element(coeffs, implicit)
    explicit = explicit_basis(a)
    return SA.coeffs(SA.coeffs(a, explicit), explicit, target)
end
