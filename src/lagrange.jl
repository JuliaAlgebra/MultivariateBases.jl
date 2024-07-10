struct ImplicitLagrangeBasis
end

struct LagrangePolynomial{T,V<:AbstractVector{T}}
    point::V
end

struct LagrangeBasis{T,P<:AbstractVector{T},V<:AbstractVector{P}}
    points::V
end
