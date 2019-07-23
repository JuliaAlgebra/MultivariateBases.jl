using Test

using MultivariateBases

@testset "Monomial" begin
    include("monomial.jl")
end
@testset "Fixed" begin
    include("fixed.jl")
end
@testset "Scaled" begin
    include("scaled.jl")
end
