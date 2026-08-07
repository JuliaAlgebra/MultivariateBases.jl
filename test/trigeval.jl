module TestTrigEval

using Test
using LinearAlgebra
using DynamicPolynomials
import StarAlgebras as SA
import MultivariateBases as MB
import TrigPolys

# A `2n+1`-point grid that matches `TrigPolys.evaluate`'s canonical layout:
# `θ_k = 2π k / (2n+1)` for `k = 0, …, 2n`. The univariate `Trigonometric`
# basis is parameterised by `value = cos θ`, so we wrap each `cos θ_k` in a
# 1-element vector to mimic how `LagrangeBasis` stores multivariate points.
_grid_points(n::Int) = [Float64[cos(2π * k / (2n + 1))] for k in 0:(2n)]

_grid_thetas(n::Int) = [2π * k / (2n + 1) for k in 0:(2n)]

# Hand-built reference: the trig polynomial
# `p(θ) = a₀ + Σₖ (c_k cos(kθ) + s_k sin(kθ))` evaluated at each grid `θ_k`.
# Coefficients are in `MultivariateBases`' interleaved order
# `[a₀, c₁, s₁, c₂, s₂, …, c_n, s_n]`.
function _eval_reference(coefs::AbstractVector, thetas::AbstractVector)
    n_coef = length(coefs)
    @assert isodd(n_coef)
    d = div(n_coef - 1, 2)
    out = zeros(eltype(coefs), length(thetas))
    for (i, θ) in pairs(thetas)
        out[i] = coefs[1]
        for k in 1:d
            out[i] += coefs[2k] * cos(k * θ) + coefs[2k+1] * sin(k * θ)
        end
    end
    return out
end

function test_construction_and_size()
    n = 2
    n_coef = 2n + 1
    pts = _grid_points(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    @test M isa MB.TrigEvalMatrix{Float64}
    @test eltype(M) === Float64
    @test size(M) == (length(pts), n_coef)
    @test size(M, 1) == length(pts)
    @test size(M, 2) == n_coef
end

function test_dense_materialization()
    # The `.dense` field is built at construction; `M[i, j]` and `Matrix(M)`
    # must forward to it, and a manual recurrence on the points must give
    # the same answer.
    n = 3
    n_coef = 2n + 1
    pts = _grid_points(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    den = Matrix(M)
    @test den == M.dense
    @test size(den) == size(M)
    for i in 1:size(M, 1), j in 1:size(M, 2)
        @test M[i, j] == den[i, j]
    end
    # Constant column is 1; first non-constant column equals each point value
    # (degree-1 univariate polynomial is `value`).
    @test all(M[:, 1] .== 1.0)
    @test M[:, 2] == [only(p) for p in pts]
end

function test_index_style_and_strides()
    # Strided forwarding is required so `view(M, j, :)` plus
    # `LinearAlgebra.dot` dispatches to BLAS rather than the generic loop.
    n = 4
    n_coef = 2n + 1
    pts = _grid_points(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    @test IndexStyle(M) == IndexLinear()
    @test IndexStyle(typeof(M)) == IndexLinear()
    @test strides(M) == strides(M.dense)
    @test Base.elsize(typeof(M)) == sizeof(Float64)
    # `unsafe_convert` must point at the same memory as the underlying dense.
    @test Base.unsafe_convert(Ptr{Float64}, M) ===
          Base.unsafe_convert(Ptr{Float64}, M.dense)
    # Linear indexing reaches the same elements as cartesian.
    @test M[1] == M[1, 1]
    @test M[length(pts)+1] == M[1, 2]
end

function test_fft_mul_matches_reference()
    # The FFT path (`M * x` for `x::AbstractVector{<:Real}`) routes through
    # `TrigPolys.evaluate` after reordering MB's interleaved layout into
    # TrigPolys' separated layout. On the canonical grid it must match the
    # hand-evaluated trig polynomial.
    n = 5
    n_coef = 2n + 1
    pts = _grid_points(n)
    thetas = _grid_thetas(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    coefs = [0.7, -0.4, 0.2, 1.1, -0.3, 0.05, 0.9, -1.2, 0.6, 0.1, -0.5]
    @assert length(coefs) == n_coef
    expected = _eval_reference(coefs, thetas)
    y_alloc = M * coefs
    @test y_alloc ≈ expected
    y_inplace = similar(y_alloc)
    LinearAlgebra.mul!(y_inplace, M, coefs)
    @test y_inplace ≈ expected
end

function test_fft_mul_matrix_rhs()
    # Matrix-RHS `mul!` iterates columns; each column must match `M * col`.
    n = 4
    n_coef = 2n + 1
    pts = _grid_points(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    k = 3
    X = randn(n_coef, k)
    Y = Matrix{Float64}(undef, size(M, 1), k)
    LinearAlgebra.mul!(Y, M, X)
    for j in 1:k
        @test Y[:, j] ≈ M * X[:, j]
    end
end

function test_adjoint_mul_is_real_adjoint()
    # `mul!(x, M', y)` routes through `TrigPolys.evaluateT`. We don't have a
    # closed form for the unpadded adjoint, so verify the inner-product
    # identity `⟨Mx, y⟩ == ⟨x, M' y⟩` directly.
    n = 4
    n_coef = 2n + 1
    pts = _grid_points(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    x = randn(n_coef)
    y = randn(size(M, 1))
    Mx = M * x
    Mty = M' * y
    @test Mty isa Vector{Float64}
    @test length(Mty) == n_coef
    @test dot(Mx, y) ≈ dot(x, Mty)
    # In-place adjoint
    Mty_ip = similar(Mty)
    LinearAlgebra.mul!(Mty_ip, M', y)
    @test Mty_ip ≈ Mty
end

# A minimal symbolic-coefficient type. Mirrors what MOI's
# `ScalarAffineFunction` looks like to `Matrix(M) * x`: an opaque element
# closed under `Float64 * _` and `_ + _`. The FFT path is gated on
# `<:Real`, so this exercises the dense fallback.
struct Symb
    label::Symbol
end
Base.zero(::Type{Symb}) = Symb(:zero)
Base.zero(x::Symb) = zero(Symb)
Base.iszero(::Symb) = false
Base.:+(::Symb, ::Symb) = Symb(:sum)
Base.:*(::Float64, ::Symb) = Symb(:scaled)
Base.:*(::Symb, ::Float64) = Symb(:scaled)
Base.promote_rule(::Type{Symb}, ::Type{Float64}) = Symb

function test_symbolic_fallback_uses_dense()
    # When the RHS isn't `<:Real`, the generic `*` overload should
    # materialize the dense matrix and hand off to `Matrix * x`. This is
    # what the SumOfSquares SOS bridge relies on when the constraint
    # function coefficients are `MOI.ScalarAffineFunction`.
    n = 2
    n_coef = 2n + 1
    pts = _grid_points(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    x = Symb[Symb(Symbol("c$i")) for i in 1:n_coef]
    y = M * x
    @test y isa AbstractVector{Symb}
    @test length(y) == size(M, 1)
    # Each output is a sum (across the row) of scaled symbolic terms.
    @test all(yi -> yi isa Symb, y)
end

function test_transformation_to_dispatch()
    # The new `transformation_to(::SubBasis{Trigonometric}, ::LagrangeBasis)`
    # method should fire and return a `TrigEvalMatrix` — NOT the generic
    # dense `Matrix{T}` built by the fallback at `lagrange.jl:119`.
    @polyvar t
    n = 3
    monos = monomials(t, 0:(2n))
    sub = MB.SubBasis{MB.Trigonometric}(monos)
    pts = _grid_points(n)
    lag = MB.LagrangeBasis([t], pts)
    U = MB.transformation_to(sub, lag)
    @test U isa MB.TrigEvalMatrix
    @test size(U) == (length(lag), length(sub))
end

function test_transformation_to_promote_op()
    # `Base.promote_op(MB.transformation_to, G, B)` must resolve the matrix
    # type at compile time so SumOfSquares' `LowRankBridge`'s
    # `added_constrained_variable_types` can declare a concrete factor type.
    @polyvar t
    n = 2
    monos = monomials(t, 0:(2n))
    sub = MB.SubBasis{MB.Trigonometric}(monos)
    pts = _grid_points(n)
    lag = MB.LagrangeBasis([t], pts)
    MT = Base.promote_op(MB.transformation_to, typeof(sub), typeof(lag))
    @test MT <: MB.TrigEvalMatrix{Float64}
    @test MT !== Any
    # Same query for a non-Trig basis should drop to the generic
    # `Matrix{Float64}` builder.
    sub_mono = MB.SubBasis{MB.Monomial}(monos)
    MT_mono =
        Base.promote_op(MB.transformation_to, typeof(sub_mono), typeof(lag))
    @test MT_mono <: AbstractMatrix{Float64}
    @test !(MT_mono <: MB.TrigEvalMatrix)
end

function test_row_view_is_strided_blas_path()
    # `view(M, j, :)` should be a strided `SubArray`. The whole point of
    # `IndexStyle = IndexLinear` + `strides` + `unsafe_convert` is that
    # downstream `dot(::SubArray, ::AbstractVector)` reaches BLAS rather than
    # falling back to the generic Julia loop.
    n = 4
    n_coef = 2n + 1
    pts = _grid_points(n)
    M = MB.TrigEvalMatrix{Float64}(pts, n_coef)
    row = view(M, 2, :)
    @test row isa SubArray
    @test parent(row) === M
    @test Vector(row) == M.dense[2, :]
    # Inner product equals the dense reference.
    v = randn(n_coef)
    @test dot(row, v) ≈ dot(M.dense[2, :], v)
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestTrigEval.runtests()
