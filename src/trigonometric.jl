"""
    struct Trigonometric <: AbstractMonomialIndexed end

Univariate trigonometric basis is
```
a0 + a1 cos(ωt) + a2 sin(ωt) + a3 cos(2ωt) + a4 sin(2ωt)
```
"""
struct Trigonometric <: AbstractMultipleOrthogonal end

_cos_id(d) = iszero(d) ? 0 : 2d - 1
_sin_id(d) = 2d
_is_cos(d) = isodd(d)
_is_sin(d) = d > 0 && iseven(d)
_id(d) = div(d + 1, 2)

# https://en.wikipedia.org/wiki/Chebyshev_polynomials#Properties
# Using
# sin(a + b) = sin(a) cos(b) + cos(a) sin(b)
# sin(a - b) = sin(a) cos(b) - cos(a) sin(b)
# If a > b
# sin(a) cos(b) = sin(a + b) + sin(a - b)
# sin(a) cos(b) = sin(a + b) + sin(a - b)
function univariate_mul!(::Type{Trigonometric}, exps, coefs, var, a, b)
    @assert !iszero(a)
    @assert !iszero(b)
    da = _id(a)
    db = _id(b)
    for i in eachindex(exps)
        if _is_cos(a) == _is_cos(b)
            # Chebyshev first kind
            push!(exps, _increment(exps[i], _cos_id(da + db), var))
            exps[i] = _increment!(exps[i], _cos_id(abs(da - db)), var)
            coefs[i] = MA.operate!!(/, coefs[i], 2)
            push!(coefs, MA.copy_if_mutable(coefs[i]))
            # cos(a + b) = cos(a) cos(b) - sin(a) sin(b)
            # cos(a - b) = cos(a) cos(b) + sin(a) sin(b)
            # cos(a - b) - cos(a + b)
            if _is_sin(a)
                coefs[end] = MA.operate!!(*, coefs[end], -1)
            end
        else
            if _is_cos(a)
                da, db = db, da
            end
            # sin(da) * cos(db)
            if da == db
                # sin(da) * cos(da) = sin(2da) / 2
                exps[i] = _increment!(exps[i], _cos_id(da + db), var)
                coefs[i] = MA.operate!!(/, coefs[i], 2)
            else
                # Using
                # sin(a + b) = sin(a) cos(b) + cos(a) sin(b)
                # sin(a - b) = sin(a) cos(b) - cos(a) sin(b)
                # If a > b
                # sin(a) cos(b) = (sin(a + b) + sin(a - b)) / 2
                # If a < b
                # sin(a) cos(b) = (sin(b + a) - sin(b - a)) / 2
                push!(exps, _increment(exps[i], _sin_id(da + db), var))
                exps[i] = _increment!(exps[i], _sin_id(abs(da - db)), var)
                coefs[i] = MA.operate!!(/, coefs[i], 2)
                push!(coefs, MA.copy_if_mutable(coefs[i]))
                if da < db
                    coefs[i] = MA.operate!!(*, coefs[i], -1)
                end
            end
        end
    end
    return
end

function degree_one_univariate_polynomial(::Type{Trigonometric}, variable)
    MA.@rewrite(variable + 0)
end

function recurrence_eval(
    ::Type{Trigonometric},
    previous::AbstractVector,
    value,
    degree,
)
    d = _id(degree)
    if _is_cos(degree)
        # Chebyshev first order
        return 2 * value * previous[_cos_id(d-1)+1] - previous[_cos_id(d-2)+1]
    else
        return sqrt(1 - previous[degree]^2)
    end
end

function _promote_coef(::Type{T}, ::Type{Trigonometric}) where {T}
    return _promote_div(T)
end

# FIXME The cos part is, like Chebysev, maybe the sin part too ? We should do better here, this is just a stopgap
even_odd_separated(::Type{Trigonometric}) = false

import TrigPolys

"""
    TrigEvalMatrix{T,P<:AbstractVector} <: AbstractMatrix{T}

Change-of-basis matrix from a univariate `SubBasis{Trigonometric}` (rows
indexing samples / columns indexing the canonical trig basis
`[1, cos(ω), sin(ω), cos(2ω), sin(2ω), ...]` in `MultivariateBases`'
interleaved order) to a `LagrangeBasis` whose nodes are `points`.

`mul!`/`Base.*` use `TrigPolys.evaluate` / `TrigPolys.evaluateT`, i.e. FFTs
on the canonical TrigPolys uniform grid.

TODO: this currently ignores `points` for the fast path — we assume they
are the `2n+1` equispaced grid points TrigPolys uses internally. Calling
`mul!` on points that do not match that grid will give silently wrong
results. A proper implementation should detect mismatch and fall back to
the dense Vandermonde matrix built via `recurrence_eval`.
"""
struct TrigEvalMatrix{T,P<:AbstractVector} <: AbstractMatrix{T}
    points::P
    n_coef::Int
    # Dense materialization computed at construction. Building it costs
    # `O(n d)` (Chebyshev-style recurrence per row), the same total work as
    # one BLAS matmul would do on a stored matrix. Holding it makes per-element
    # `getindex` `O(1)` — required so the downstream `LRO.Factorization`'s
    # `dot` chain (which iterates `factor[k]` for `factor::SubArray{T,1,<:TrigEvalMatrix}`)
    # stays cheap. `mul!` on `AbstractVector{<:Real}` / `AbstractMatrix` columns
    # still routes through `TrigPolys.evaluate`/`evaluateT`, so when downstream
    # code dispatches at the parent-matrix level (batched FFT in BM, see TODO
    # above) the FFT path remains live.
    dense::Matrix{T}
end

function _materialize_trig_eval(
    ::Type{T},
    points::AbstractVector,
    n_coef::Integer,
) where {T}
    n_pts = length(points)
    out = Matrix{T}(undef, n_pts, n_coef)
    for i in 1:n_pts
        # `LagrangeBasis` stores each point as a 1-element `AbstractVector`
        # (univariate); the basis-recurrence expects the scalar value.
        val = T(only(points[i]))
        if n_coef >= 1
            out[i, 1] = one(T)
        end
        if n_coef >= 2
            out[i, 2] = degree_one_univariate_polynomial(Trigonometric, val)
        end
        for d in 2:(n_coef-1)
            out[i, d+1] =
                recurrence_eval(Trigonometric, view(out, i, 1:d), val, d)
        end
    end
    return out
end

function TrigEvalMatrix{T}(
    points::P,
    n_coef::Integer,
) where {T,P<:AbstractVector}
    return TrigEvalMatrix{T,P}(
        points,
        n_coef,
        _materialize_trig_eval(T, points, n_coef),
    )
end

Base.size(M::TrigEvalMatrix) = (length(M.points), M.n_coef)

# Forward the strided-array interface to the underlying dense storage so that
# `view(M, j, :)` is a strided `SubArray` that BLAS' `gemv!` accepts. Without
# this, downstream `LinearAlgebra.dot(::Factorization{T,<:SubArray{T,1,<:TrigEvalMatrix}}, ...)`
# in `LowRankOpt.BurerMonteiro` would fall back to the generic Julia loop with
# allocations per dot — orders of magnitude slower than BLAS.
Base.IndexStyle(::Type{<:TrigEvalMatrix}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(M::TrigEvalMatrix, i::Int) = M.dense[i]
Base.@propagate_inbounds Base.getindex(
    M::TrigEvalMatrix,
    i::Integer,
    j::Integer,
) = M.dense[i, j]
Base.strides(M::TrigEvalMatrix) = strides(M.dense)
function Base.unsafe_convert(::Type{Ptr{T}}, M::TrigEvalMatrix{T}) where {T}
    return Base.unsafe_convert(Ptr{T}, M.dense)
end
Base.elsize(::Type{<:TrigEvalMatrix{T}}) where {T} = sizeof(T)

# Reorder a coefficient vector laid out in MultivariateBases' interleaved order
# `[a0, c1, s1, c2, s2, ...]` into TrigPolys' layout `[a0, c1, ..., cn, s1, ..., sn]`.
function _trig_coefs_to_trigpoly!(out::AbstractVector, x::AbstractVector)
    n_x = length(x)
    n_out = length(out)
    @assert isodd(n_x)
    @assert isodd(n_out)
    d = div(n_x - 1, 2)
    n = div(n_out - 1, 2)
    @assert n >= d
    out[1] = x[1]
    @inbounds for k in 1:d
        out[1+k] = x[2k]
        out[1+n+k] = x[2k+1]
    end
    return out
end

function _trigpoly_to_trig_coefs!(out::AbstractVector, a::AbstractVector)
    n_out = length(out)
    n_a = length(a)
    @assert isodd(n_out)
    @assert isodd(n_a)
    d = div(n_out - 1, 2)
    n = div(n_a - 1, 2)
    @assert n >= d
    out[1] = a[1]
    @inbounds for k in 1:d
        out[2k] = a[1+k]
        out[2k+1] = a[1+n+k]
    end
    return out
end

function _trigpoly_for_eval(M::TrigEvalMatrix{T}, x::AbstractVector) where {T}
    n_coef = M.n_coef
    @assert length(x) == n_coef
    @assert isodd(n_coef)
    grid_d = div(length(M.points) - 1, 2)
    n = max(div(n_coef - 1, 2), grid_d)
    a = zeros(T, 2n + 1)
    # `a` is laid out as `[a0, c1..cn, s1..sn]`; padding for the high
    # frequencies stays at zero (mirroring `vectorized_pad_to`).
    _trig_coefs_to_trigpoly!(a, x)
    return TrigPolys.TrigPoly(a)
end

function LinearAlgebra.mul!(
    y::AbstractVector{<:Real},
    M::TrigEvalMatrix{T},
    x::AbstractVector{<:Real},
) where {T}
    p = _trigpoly_for_eval(M, x)
    copyto!(y, TrigPolys.evaluate(p))
    return y
end

function Base.:*(M::TrigEvalMatrix{T}, x::AbstractVector{<:Real}) where {T}
    y = Vector{T}(undef, size(M, 1))
    return LinearAlgebra.mul!(y, M, x)
end

# Fallback when `x` carries symbolic / MOI-style entries (e.g. during
# `SA.coeffs(::AlgebraElement{<:MOI.ScalarAffineFunction}, ...)` in the SOS
# constraint bridge). The FFT pad/unpad path needs numeric storage, so we
# materialize a dense matrix and hand off to the generic `*`.
function Base.:*(M::TrigEvalMatrix, x::AbstractVector)
    return Matrix(M) * x
end

function LinearAlgebra.mul!(
    Y::AbstractMatrix,
    M::TrigEvalMatrix{T},
    X::AbstractMatrix,
) where {T}
    @assert size(Y, 1) == size(M, 1)
    @assert size(M, 2) == size(X, 1)
    for j in axes(X, 2)
        @views LinearAlgebra.mul!(Y[:, j], M, X[:, j])
    end
    return Y
end

# Adjoint operator: `evaluateT` is the adjoint of `evaluate`.
function LinearAlgebra.mul!(
    x::AbstractVector,
    Madj::LinearAlgebra.Adjoint{<:Any,<:TrigEvalMatrix{T}},
    y::AbstractVector,
) where {T}
    M = parent(Madj)
    @assert length(x) == M.n_coef
    @assert length(y) == size(M, 1)
    # `a` is laid out as `[a0, c1..cn, s1..sn]`; `_trigpoly_to_trig_coefs!`
    # unpads it into MB's interleaved layout.
    _trigpoly_to_trig_coefs!(x, TrigPolys.evaluateT(y))
    return x
end

function Base.:*(
    Madj::LinearAlgebra.Adjoint{<:Any,<:TrigEvalMatrix{T}},
    y::AbstractVector,
) where {T}
    x = Vector{T}(undef, size(Madj, 1))
    return LinearAlgebra.mul!(x, Madj, y)
end

# Override `transformation_to(::SubBasis{Trigonometric}, ::LagrangeBasis)` so
# that the bridge layer carries an FFT-backed `AbstractMatrix` instead of the
# dense Vandermonde matrix built by the generic `transformation_to` in
# `lagrange.jl`. The new lazy matrix multiplies via TrigPolys' FFT, so the
# downstream `LowRankOpt.BurerMonteiro` `mul!` calls run in `O(d log d)`.
