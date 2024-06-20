struct MultiBasis{T,I,B<:SA.ExplicitBasis{T,I}} <: SA.ExplicitBasis{T,I}
    bases::Vector{B}
end

Base.length(basis::MultiBasis) = length(first(basis.bases))
