module MultivariateBasesMathOptInterfaceExt

import MathOptInterface as MOI
import MultivariateBases as MB

function MB._matrix_eltype(
    ::Type{F},
) where {T,F<:Union{MOI.ScalarAffineFunction{T},MOI.ScalarQuadraticFunction{T}}}
    return T
end

end
