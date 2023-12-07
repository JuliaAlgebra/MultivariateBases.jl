abstract type MonomialIndexedBasis{T,M<:MP.AbstractMonomial} <: SA.ImplicitBasis{T,M} end

export ChebyshevFirstKind, ChebyshevFirstKindBasis

struct ChebyshevFirstKind{M<:MP.AbstractMonomial}
    monomial::M
end
ChebyshevFirstKind(v::MP.AbstractVariable) = ChebyshevFirstKind(MP.monomial(v))
function _show(io::IO, mime::MIME, cheby::ChebyshevFirstKind)
    print(io, "Cheby1st(")
    show(io, mime, cheby.monomial)
    print(io, ")")
end
function Base.show(io::IO, mime::MIME"text/plain", cheby::ChebyshevFirstKind)
    _show(io, mime, cheby)
end
function Base.show(io::IO, cheby::ChebyshevFirstKind)
    show(io, MIME"text/plain"(), cheby)
end

struct ChebyshevFirstKindBasis{M} <: MonomialIndexedBasis{ChebyshevFirstKind{M},M} end
function Base.getindex(::ChebyshevFirstKindBasis{M}, mono::M) where {M}
    return ChebyshevFirstKind(mono)
end

function Base.zero(::Type{ChebyshevFirstKind{M}}) where {M}
    basis = ChebyshevFirstKindBasis{M}()
    return SA.AlgebraElement(
        SA.SparseCoefficients(M[], Rational{Int}[]),
        SA.StarAlgebra(
            ChebyshevFirstKind(MP.constant_monomial(M)),
            basis,
            SA.LazyMStructure(basis),
        ),
    )
end

Base.zero(a::ChebyshevFirstKind) = zero(typeof(a))

function Base.:*(a::ChebyshevFirstKind{M}, b::ChebyshevFirstKind{M}) where {M}
    return MA.operate!(MA.add_mul, zero(a), 1 // 1, a, b)
end

function MA.operate!(
    ::typeof(MA.add_mul),
    result::SA.AlgebraElement{A,V,<:SA.SparseCoefficients},
    α::V,
    a::ChebyshevFirstKind,
    b::ChebyshevFirstKind,
) where {A,V}
    # https://en.wikipedia.org/wiki/Chebyshev_polynomials#Properties
    # T_n * T_m = T_{n + m} / 2 + T_{|n - m|} / 2
    coeffs = result.coeffs
    push!(coeffs.basis_elements, MP.constant_monomial(a.monomial * b.monomial))
    push!(coeffs.values, α)
    i = length(coeffs.values)
    vars_b = MP.variables(b.monomial)
    var_state_b = iterate(vars_b)
    for var_a in MP.variables(a.monomial)
        if isnothing(var_state_b)
            break
        end
        var_b, state_b = var_state_b
        if var_b > var_a
            var_state_b = iterate(vars_b, state_b)
            for j in i:length(coeffs.basis_elements)
                coeffs.basis_elements[j] = MA.mul!!(coeffs.basis_elements[j], var_b^MP.degree(b.monomial, var_b))
            end
        elseif var_a > var_b
            for j in i:length(coeffs.basis_elements)
                coeffs.basis_elements[j] = MA.mul!!(coeffs.basis_elements[j], var_a^MP.degree(a.monomial, var_a))
            end
        else
            n = length(coeffs.basis_elements)
            d_a = MP.degree(a.monomial, var_a)
            d_b = MP.degree(b.monomial, var_b)
            for j in i:n
                push!(coeffs.basis_elements, coeffs.basis_elements[j] * var_a^abs(d_a - d_b))
                coeffs.basis_elements[j] = MA.mul!!(coeffs.basis_elements[j], var_a^(d_a + d_b))
                coeffs.values[j] = MA.operate!!(/, coeffs.values[j], 2)
                push!(coeffs.values, MA.copy_if_mutable(coeffs.values[j]))
            end
        end
    end
    return result
end

abstract type AbstractChebyshevBasis{P} <: AbstractGegenbauerBasis{P} end

function MP.polynomial_type(::Type{<:AbstractChebyshevBasis}, V::Type)
    return MP.polynomial_type(V, Float64)
end

reccurence_first_coef(::Type{<:AbstractChebyshevBasis}, degree) = 2
reccurence_third_coef(::Type{<:AbstractChebyshevBasis}, degree) = -1
reccurence_deno_coef(::Type{<:AbstractChebyshevBasis}, degree) = 1

"""
    struct ChebyshevBasisFirstKind{P} <: AbstractChebyshevBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\frac{1}{\\sqrt{1 - x^2}}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevBasisFirstKind{P} <: AbstractChebyshevBasis{P}
    polynomials::Vector{P}
end

const ChebyshevBasis{P} = ChebyshevBasisFirstKind{P}

function degree_one_univariate_polynomial(
    ::Type{<:ChebyshevBasisFirstKind},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(variable + 0)
end

function _scalar_product_function(::Type{<:ChebyshevBasisFirstKind}, i::Int)
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
    struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
        polynomials::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\sqrt{1 - x^2}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
    polynomials::Vector{P}
end

function degree_one_univariate_polynomial(
    ::Type{<:ChebyshevBasisSecondKind},
    variable::MP.AbstractVariable,
)
    MA.@rewrite(2variable + 0)
end

function _scalar_product_function(::Type{<:ChebyshevBasisSecondKind}, i::Int)
    if i == 0
        return π / 2
    elseif isodd(i)
        return 0
    else
        n = div(i, 2)
        return π / (2^(i + 1)) * prod(n+2:i) / factorial(n)
    end
end
