const _APL = MP.AbstractPolynomialLike
# We don't define it for all `AlgebraElement` as this would be type piracy
const _AE = SA.AlgebraElement{<:Algebra}

struct SubMStruct{B} <: SA.MultiplicativeStructure
    bases::B
end

function (m::SubMStruct)(args...)
    return SA.coeffs(
        *(getindex.(m.bases[2:end], args)...),
        m.bases[1],
    )
end

function SA.mstructure(args::Algebra...)
    if SA.basis(args[1]) isa FullBasis &&
        all(isequal(args[1]), args)
        return SA.mstructure(SA.basis(args[1]))
    end
    return SubMStruct(SA.basis.(args))
end

SA.mstructure(::FullBasis{B}) where {B} = Mul{B}()

function _preallocate_output(op, args::Vararg{Any,N}) where {N}
    T = MA.promote_operation(op, SA._coeff_type.(args)...)
    return zero(T, algebra(implicit_basis(SA.basis(args[1]))))
end

for op in [:+, :-, :*]
    @eval begin
        function MA.promote_operation(
            ::typeof($op),
            ::Type{P},
            ::Type{Q},
        ) where {P<:_APL,Q<:_AE}
            return MA.promote_operation($op, P, MP.polynomial_type(Q))
        end
        Base.$op(p::_APL, q::_AE) = $op(p, MP.polynomial(q))
        function MA.promote_operation(
            ::typeof($op),
            ::Type{P},
            ::Type{Q},
        ) where {P<:_AE,Q<:_APL}
            return MA.promote_operation($op, MP.polynomial_type(P), Q)
        end
        Base.$op(p::_AE, q::_APL) = $op(MP.polynomial(p), q)
        # Break ambiguity between the two defined below and the generic one in SA
        function MA.promote_operation(
            ::typeof($op),
            ::Type{P},
            ::Type{Q},
        ) where {P<:_AE,Q<:_AE}
            return SA.algebra_promote_operation($op, P, Q)
        end
        function Base.$op(p::_AE, q::_AE)
            return MA.operate_to!(_preallocate_output($op, p, q), $op, p, q)
        end
    end
end
for op in [:+, :-]
    @eval begin
        function MA.promote_operation(
            ::typeof($op),
            ::Type{P},
            ::Type{Q},
        ) where {P,Q<:_AE}
            I = MA.promote_operation(implicit, Q)
            return MA.promote_operation(
                $op,
                constant_algebra_element_type(
                    MA.promote_operation(SA.basis, I),
                    P,
                ),
                I,
            )
        end
        function Base.$op(p, q::_AE)
            i = implicit(q)
            return $op(constant_algebra_element(typeof(SA.basis(i)), p), i)
        end
        function MA.promote_operation(
            ::typeof($op),
            ::Type{P},
            ::Type{Q},
        ) where {P<:_AE,Q}
            I = MA.promote_operation(implicit, P)
            return MA.promote_operation(
                $op,
                I,
                constant_algebra_element_type(
                    MA.promote_operation(SA.basis, I),
                    Q,
                ),
            )
        end
        function Base.$op(p::_AE, q)
            i = implicit(p)
            return $op(i, constant_algebra_element(typeof(SA.basis(i)), q))
        end
    end
end

function term_element(α, p::Polynomial{B,M}) where {B,M}
    return algebra_element(
        sparse_coefficients(MP.term(α, p.monomial)),
        FullBasis{B,M}(),
    )
end

# Needed by `SymbolicWedderburn` which multiplies elements of the basis by `Int`
# We'll see if `::Number` is too restrictive
# Should be able to remove once https://github.com/kalmarek/SymbolicWedderburn.jl/issues/88 is closed
Base.:*(α::Number, p::Polynomial) = term_element(α, p)

function MA.operate!(op::Union{typeof(+),typeof(-),typeof(*)}, p::_APL, q::_AE)
    return MA.operate!(op, p, MP.polynomial(q))
end

function MA.operate_to!(
    res::MP.AbstractPolynomial,
    op::typeof(*),
    p::_AE,
    q::_APL,
)
    return MA.operate_to!(res, op, MP.polynomial(p), q)
end

function MA.operate_to!(
    res::MP.AbstractPolynomial,
    op::typeof(*),
    p::_APL,
    q::_AE,
)
    return MA.operate_to!(res, op, p, MP.polynomial(q))
end

# These are not implemented yet for arbitrary bases so we
# fall back to polynomials

function MP.substitute(
    s::MP.AbstractSubstitutionType,
    p::_AE,
    args::MP.Substitutions,
)
    return MP.substitute(s, MP.polynomial(p), args)
end

function MP.subs(p::_AE, args::MP.AbstractSubstitution...)
    return MP.substitute(MP.Subs(), p, args)
end

function (p::_AE)(args::MP.AbstractSubstitution...)
    return MP.substitute(MP.Eval(), p, args)
end

function (p::_AE)(x::NTuple{N,<:Number}) where {N}
    return (MP.polynomial(p))(x)
end

function (p::_AE)(x::AbstractVector{<:Number})
    return (MP.polynomial(p))(x)
end

(p::_AE)(x::Number...) = (MP.polynomial(p))(x...)

function MP.differentiate(p::_AE, x)
    return MP.differentiate(MP.polynomial(p), x)
end
