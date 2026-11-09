module MultivariateBasesJuMPExt

import JuMP
import MultivariateBases as MB

function MB._matrix_eltype(::Type{F}) where {F<:JuMP.AbstractJuMPScalar}
    return JuMP.value_type(F)
end

end
