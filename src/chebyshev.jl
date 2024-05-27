abstract type AbstractChebyshev <: AbstractGegenbauer end

_promote_div(::Type{I}) where {I<:Integer} = Rational{I}
_promote_div(::Type{T}) where {T<:Number} = MA.promote_operation(/, T, Int)
# Could be for instance `MathOptInterface.ScalarAffineFunction{Float64}`
# which does not support division with `Int`
_promote_div(::Type{F}) where {F} = F

function MP.polynomial_type(
    ::Type{Polynomial{B,M}},
    ::Type{T},
) where {B<:AbstractChebyshev,M,T}
    return MP.polynomial_type(M, _promote_div(T))
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
# T_n * T_m = T_{n + m} / 2 + T_{|n - m|} / 2
function (::Mul{Chebyshev})(a::MP.AbstractMonomial, b::MP.AbstractMonomial)
    terms = [MP.term(1 // 1, MP.constant_monomial(a * b))]
    vars_b = MP.variables(b)
    var_state_b = iterate(vars_b)
    for var_a in MP.variables(a)
        if isnothing(var_state_b)
            break
        end
        var_b, state_b = var_state_b
        if var_b > var_a
            var_state_b = iterate(vars_b, state_b)
            for i in eachindex(terms)
                terms[i] = MA.mul!!(terms[i], var_b^MP.degree(b, var_b))
            end
        elseif var_a > var_b
            for i in eachindex(terms)
                terms[i] = MA.mul!!(terms[i], var_a^MP.degree(a, var_a))
            end
        else
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
        end
    end
    return MP.polynomial!(terms)
end

function SA.coeffs(coeffs, basis::FullBasis{Chebyshev}, ::FullBasis{Monomial})
    res = zero(MP.polynomial_type(typeof(basis), valtype(coeffs)))
    for (k, v) in SA.nonzero_pairs(coeffs)
        MA.operate!(SA.UnsafeAddMul(*), res, v, MP.polynomial(basis[k]))
    end
    MA.operate!(SA.canonical, res)
    return res
end

function degree_one_univariate_polynomial(
    ::Type{Chebyshev},
    variable::MP.AbstractVariable,
)
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

function degree_one_univariate_polynomial(
    ::Type{ChebyshevSecondKind},
    variable::MP.AbstractVariable,
)
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
