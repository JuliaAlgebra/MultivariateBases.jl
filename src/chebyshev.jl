abstract type AbstractChebyshev <: AbstractGegenbauer end

_promote_div(::Type{I}) where {I<:Integer} = Rational{I}
_promote_div(::Type{T}) where {T<:Number} = MA.promote_operation(/, T, Int)
# Could be for instance `MathOptInterface.ScalarAffineFunction{Float64}`
# which does not support division with `Int`
_promote_div(::Type{F}) where {F} = _float(F)

function _promote_coef(::Type{T}, ::Type{<:AbstractChebyshev}) where {T}
    return _promote_div(T)
end

reccurence_first_coef(::Type{<:AbstractChebyshev}, degree) = 2
reccurence_third_coef(::Type{<:AbstractChebyshev}, degree) = -1
reccurence_deno_coef(::Type{<:AbstractChebyshev}, degree) = 1

"""
    struct ChebyshevFirstKind <: AbstractChebyshev end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\frac{1}{\\sqrt{1 - x^2}}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevFirstKind <: AbstractChebyshev end
const Chebyshev = ChebyshevFirstKind

# https://en.wikipedia.org/wiki/Chebyshev_polynomials#Properties
# `T_n * T_m = T_{n + m} / 2 + T_{|n - m|} / 2`
function (::Mul{Chebyshev})(a::MP.AbstractMonomial, b::MP.AbstractMonomial)
    terms = [MP.term(1 // 1, MP.constant_monomial(a * b))]
    vars_a = MP.variables(a)
    var_state_a = iterate(vars_a)
    vars_b = MP.variables(b)
    var_state_b = iterate(vars_b)
    while !isnothing(var_state_a) || !isnothing(var_state_b)
        if isnothing(var_state_a) ||
           (!isnothing(var_state_b) && var_state_b[1] > var_state_a[1])
            var_b, state_b = var_state_b
            for i in eachindex(terms)
                terms[i] = MA.mul!!(terms[i], var_b^MP.degree(b, var_b))
            end
            var_state_b = iterate(vars_b, state_b)
        elseif isnothing(var_state_b) ||
               (!isnothing(var_state_a) && var_state_a[1] > var_state_b[1])
            var_a, state_a = var_state_a
            for i in eachindex(terms)
                terms[i] = MA.mul!!(terms[i], var_a^MP.degree(a, var_a))
            end
            var_state_a = iterate(vars_a, state_a)
        else
            var_a, state_a = var_state_a
            var_b, state_b = var_state_b
            d_a = MP.degree(a, var_a)
            d_b = MP.degree(b, var_b)
            I = eachindex(terms)
            for i in I
                mono = MP.monomial(terms[i]) * var_a^(d_a + d_b)
                terms[i] = MA.mul!!(terms[i], var_a^abs(d_a - d_b))
                terms[i] = MA.operate!!(/, terms[i], 2)
                α = MA.copy_if_mutable(MP.coefficient(terms[i]))
                push!(terms, MP.term(α, mono))
            end
            var_state_a = iterate(vars_a, state_a)
            var_state_b = iterate(vars_b, state_b)
        end
    end
    return sparse_coefficients(MP.polynomial!(terms))
end

function _add_mul_scalar_vector!(res, ::SubBasis, scalar, vector)
    return MA.operate!(MA.add_mul, res, scalar, vector)
end

function _add_mul_scalar_vector!(res, ::FullBasis, scalar, vector)
    for (k, v) in SA.nonzero_pairs(vector)
        SA.unsafe_push!(res, k, scalar * v)
    end
end

function SA.coeffs!(
    res,
    cfs,
    source::MonomialIndexedBasis{<:AbstractMultipleOrthogonal},
    target::MonomialIndexedBasis{Monomial},
)
    MA.operate!(zero, res)
    for (k, v) in SA.nonzero_pairs(cfs)
        _add_mul_scalar_vector!(res, target, v, SA.coeffs(source[k], target))
    end
    MA.operate!(SA.canonical, res)
    return res
end

function transformation_to(
    source::SubBasis{Chebyshev},
    target::SubBasis{Monomial},
)
    A = zeros(Float64, length(target), length(source))
    for (i, cheby) in enumerate(source)
        A[:, i] = SA.coeffs(algebra_element(cheby), target)
    end
    return LinearAlgebra.UpperTriangular(A)
end

function SA.coeffs(
    cfs,
    source::SubBasis{Monomial},
    target::FullBasis{Chebyshev},
)
    sub = explicit_basis_covering(target, source)
    # Need to make A square so that it's UpperTriangular
    extended = SubBasis{Monomial}(sub.monomials)
    ext = SA.coeffs(algebra_element(cfs, source), extended)
    return SA.SparseCoefficients(
        sub.monomials,
        #transformation_to(sub, extended) \ ext, # Julia v1.6 converts the matrix to the eltype of the `result` which is bad for JuMP
        LinearAlgebra.ldiv!(
            zeros(_promote_coef(eltype(ext), Chebyshev), length(sub)),
            transformation_to(sub, extended),
            ext,
        ),
    )
end

function SA.coeffs(cfs, ::FullBasis{Monomial}, target::FullBasis{Chebyshev})
    return SA.coeffs(
        SA.values(cfs),
        SubBasis{Monomial}(collect(SA.keys(cfs))),
        target,
    )
end

function degree_one_univariate_polynomial(::Type{Chebyshev}, variable)
    MA.@rewrite(variable + 0)
end

function _scalar_product_function(::Type{Chebyshev}, i::Int)
    if i == 0
        return π
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return (π / 2^i) * prod(n+1:i) / factorial(n)
    end
end

"""
    struct ChebyshevSecondKind <: AbstractChebyshevBasis end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\sqrt{1 - x^2}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevSecondKind <: AbstractChebyshev end

function degree_one_univariate_polynomial(::Type{ChebyshevSecondKind}, variable)
    MA.@rewrite(2variable + 0)
end

function _scalar_product_function(::Type{<:ChebyshevSecondKind}, i::Int)
    if i == 0
        return π / 2
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return π / (2^(i + 1)) * prod(n+2:i) / factorial(n)
    end
end
