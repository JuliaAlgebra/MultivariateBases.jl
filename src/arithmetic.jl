const _APL = MP.AbstractPolynomialLike
# We don't define it for all `AlgebraElement` as this would be type piracy
const _AE = SA.AlgebraElement{<:Algebra}

for op in [:+, :-, :*]
    @eval begin
        function MA.promote_operation(::typeof($op), ::Type{P}, ::Type{Q}) where {P<:_APL,Q<:_AE}
            return MA.promote_operation($op, P, MP.polynomial_type(Q))
        end
        Base.$op(p::_APL, q::_AE) = $op(p, MP.polynomial(q))
        function MA.promote_operation(::typeof($op), ::Type{P}, ::Type{Q}) where {P<:_AE,Q<:_APL}
            return MA.promote_operation($op, MP.polynomial_type(P), Q)
        end
        Base.$op(p::_AE, q::_APL) = $op(MP.polynomial(p), q)
        # Break ambiguity between the two defined below and the generic one in SA
        function MA.promote_operation(::typeof($op), ::Type{P}, ::Type{Q}) where {P<:_AE,Q<:_AE}
            return SA.algebra_promote_operation($op, P, Q)
        end
        function Base.$op(p::_AE, q::_AE)
            return MA.operate_to!(SA._preallocate_output($op, p, q), $op, p, q)
        end
    end
end
for op in [:+, :-]
    @eval begin
        function MA.promote_operation(::typeof($op), ::Type{P}, ::Type{Q}) where {P,Q<:_AE}
            I = MA.promote_operation(implicit, Q)
            return MA.promote_operation(
                $op,
                constant_algebra_element_type(MA.promote_operation(SA.basis, I), P),
                I,
            )
        end
        function Base.$op(p, q::_AE)
            i = implicit(q)
            return $op(constant_algebra_element(typeof(SA.basis(i)), p), i)
        end
        function MA.promote_operation(::typeof($op), ::Type{P}, ::Type{Q}) where {P<:_AE,Q}
            I = MA.promote_operation(implicit, P)
            return MA.promote_operation(
                $op,
                I,
                constant_algebra_element_type(MA.promote_operation(SA.basis, I), Q),
            )
        end
        function Base.$op(p::_AE, q)
            i = implicit(p)
            return $op(i, constant_algebra_element(typeof(SA.basis(i)), q))
        end
    end
end
