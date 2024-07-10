# Inspired by `Hypatia.jl/src/PolyUtils/realinterp.jl`

export BoxSampling

abstract type AbstractNodes{T,V} end

mutable struct BoxSampling{T,V} <: AbstractNodes{T,V}
    lower::V
    upper::V
    sample_factor::Int
    function BoxSampling(lower, upper; sample_factor = 0)
        @assert length(lower) == length(upper)
        l = float.(lower)
        u = float.(upper)
        V = promote_type(typeof(l), typeof(u))
        return new{eltype(V),V}(l, u, sample_factor)
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

struct LagrangePolynomial{T,P,V}
    variables::V
    point::P
    function LagrangePolynomial(variables, point)
        return new{eltype(point),typeof(point),typeof(variables)}(
            variables,
            point,
        )
    end
end

struct ImplicitLagrangeBasis{T,P,N<:AbstractNodes{T,P},V} <:
       SA.ImplicitBasis{LagrangePolynomial{T,P,V},Pair{V,P}}
    variables::V
    nodes::AbstractNodes{T,P}
    function ImplicitLagrangeBasis(
        variables,
        nodes::AbstractNodes{T,P},
    ) where {T,P}
        return new{T,P,typeof(nodes),typeof(variables)}(variables, nodes)
    end
end

function Base.getindex(
    implicit::ImplicitLagrangeBasis{T,P,N,V},
    subs::Pair{V,P},
) where {T,P,N,V}
    if subs.first != implicit.variables
        error(
            "Variables `$(subs.first)` do not match Lagrange basis variables `$(implicit.variables)`",
        )
    end
    return LagrangePolynomial(implicit.variables, subs.second)
end

struct LagrangeBasis{T,P,U<:AbstractVector{P},V} <:
       SA.ExplicitBasis{LagrangePolynomial{T,P,V},Int}
    variables::V
    points::U
    function LagrangeBasis(variables, points::AbstractVector)
        P = eltype(points)
        return new{eltype(P),P,typeof(points),typeof(variables)}(
            variables,
            points,
        )
    end
end

Base.length(basis::LagrangeBasis) = length(basis.points)
MP.nvariables(basis::LagrangeBasis) = length(basis.variables)
MP.variables(basis::LagrangeBasis) = basis.variables

function explicit_basis_type(
    ::Type{<:ImplicitLagrangeBasis{T,_P,N,V}},
) where {T,_P,N,V}
    points = _eachcol(ones(T, 1, 1))
    P = eltype(points)
    return LagrangeBasis{eltype(P),P,typeof(points),V}
end

function eval_basis!(
    univariate_buffer,
    result,
    basis::SubBasis{B},
    variables,
    values,
) where {B}
    for v in MP.variables(basis)
        if !(v in variables)
            error(
                "Cannot evaluate `$basis` as its variable `$v` is not part of the variables `$variables` of the `LagrangeBasis`",
            )
        end
    end
    for i in eachindex(values)
        l = MP.maxdegree(basis.monomials, variables[i]) + 1
        univariate_eval!(B, view(univariate_buffer, 1:l, i), values[i])
    end
    for i in eachindex(basis)
        result[i] = one(eltype(result))
        for j in eachindex(values)
            d = MP.degree(basis.monomials[i], variables[j])
            result[i] = MA.operate!!(*, result[i], univariate_buffer[d+1, j])
        end
    end
    return result
end

function transformation_to(basis::SubBasis, lag::LagrangeBasis{T}) where {T}
    # To avoid allocating this too often, we allocate it once here
    # and reuse it for each sample
    univariate_buffer = Matrix{T}(undef, length(basis), MP.nvariables(lag))
    V = Matrix{T}(undef, length(lag), length(basis))
    for i in eachindex(lag)
        eval_basis!(
            univariate_buffer,
            view(V, i, :),
            basis,
            MP.variables(lag),
            lag.points[i],
        )
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

if VERSION >= v"1.10"
    _eachcol(x) = eachcol(x)
    _column_norm() = LinearAlgebra.ColumnNorm()
else
    # It is a `Base.Generator` so not an `AbstractVector`
    _eachcol(x) = collect(eachcol(x))
    _column_norm() = Val(true)
end

function sample(variables, s::AbstractNodes, basis::SubBasis)
    samples = sample(s, num_samples(s.sample_factor, length(basis)))
    full = LagrangeBasis(variables, _eachcol(samples))
    V = transformation_to(basis, full)
    F = LinearAlgebra.qr!(Matrix(V'), _column_norm())
    kept_indices = F.p[1:length(basis)]
    return LagrangeBasis(variables, _eachcol(samples[:, kept_indices]))
end

function explicit_basis_covering(
    implicit::ImplicitLagrangeBasis,
    basis::SubBasis,
)
    return sample(implicit.variables, implicit.nodes, basis)
end

function SA.coeffs(coeffs, source::SubBasis, target::LagrangeBasis)
    return transformation_to(source, target) * coeffs
end

function SA.coeffs(coeffs, implicit::FullBasis, target::LagrangeBasis)
    a = algebra_element(coeffs, implicit)
    explicit = explicit_basis(a)
    return SA.coeffs(SA.coeffs(a, explicit), explicit, target)
end
