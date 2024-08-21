using Test

import MutableArithmetics as MA
import StarAlgebras as SA
using MultivariateBases
const MB = MultivariateBases
using LinearAlgebra
using DynamicPolynomials

for file in readdir(@__DIR__)
    if endswith(file, ".jl") && !in(file, ["runtests.jl"])
        @testset "$file" begin
            include(joinpath(@__DIR__, file))
        end
    end
end
