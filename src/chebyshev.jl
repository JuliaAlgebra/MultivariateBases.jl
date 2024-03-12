# TODO Add to MultivariatePolynomials
Base.keytype(p::MP.AbstractPolynomial) = MP.monomial_type(p)
Base.valtype(p::MP.AbstractPolynomial) = MP.coefficient_type(p)
#Base.keys(p::MP.AbstractPolynomial) = MP.monomials(p)
Base.pairs(p::MP.AbstractPolynomial) = MP.terms(p)
function Base.similar(p::PT, ::Type{T}) where {PT<:MP._APL,T}
    return convert(MP.similar_type(PT, T), copy(p))
end
function Base.indexed_iterate(t::MP.Term, i::Int)
    if i == 1
        return MP.coefficient(t)
    else
        return MP.monomial(t)
    end
end


struct Polynomial{B,M<:MP.AbstractMonomial}
    monomial::M
    function Polynomial{B}(mono::MP.AbstractMonomial) where {B}
        return new{B,typeof(mono)}(mono)
    end
end

function Polynomial{B}(v::MP.AbstractVariable) where {B}
    return Polynomial{B}(MP.monomial(v))
end

function _algebra_element(p, ::Type{B}) where {B}
    basis = MonomialIndexedBasis{B,MP.monomial_type(p)}()
    return SA.AlgebraElement(
        p,
        SA.StarAlgebra(
            Polynomial{B}(MP.constant_monomial(p)),
            basis,
        ),
    )
end

function Base.:*(a::Polynomial{B}, b::Polynomial{B}) where {B}
    _algebra_element(Mul{B}()(a.monomial, b.monomial), B)
end

function _show(io::IO, mime::MIME, p::Polynomial{B}) where {B}
    print(io, B)
    print(io, "(")
    show(io, mime, p.monomial)
    print(io, ")")
end
function Base.show(io::IO, mime::MIME"text/plain", p::Polynomial)
    _show(io, mime, p)
end
function Base.show(io::IO, p::Polynomial)
    show(io, MIME"text/plain"(), p)
end

struct MonomialIndexedBasis{B,M<:MP.AbstractMonomial} <: SA.ImplicitBasis{Polynomial{B,M},M} end
function Base.getindex(::MonomialIndexedBasis{B,M}, mono::M) where {B,M}
    return Polynomial{B}(mono)
end
function Base.getindex(::MonomialIndexedBasis{B,M}, p::Polynomial{B,M}) where {B,M}
    return mono
end
SA.mstructure(::MonomialIndexedBasis{B}) where {B} = Mul{B}()

function Base.zero(::Type{Polynomial{B,M}}) where {B,M}
    return _algebra_element(zero(MP.polynomial_type(M, Rational{Int})), B)
end

Base.zero(p::Polynomial) = zero(typeof(p))

struct Mul{B} <: SA.MultiplicativeStructure end

export ChebyshevFirstKind, ChebyshevFirstKindBasis

struct ChebyshevFirstKind end
const Chebyshev = ChebyshevFirstKind

export Chebyshev

# https://en.wikipedia.org/wiki/Chebyshev_polynomials#Properties
# T_n * T_m = T_{n + m} / 2 + T_{|n - m|} / 2
function (::Mul{Chebyshev})(a::MP.AbstractMonomial, b::MP.AbstractMonomial)
    terms = [MP.term(1//1, MP.constant_monomial(a * b))]
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
