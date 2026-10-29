include("utils.jl")

for file in readdir(@__DIR__)
    if endswith(file, ".jl") && !in(file, ["runtests.jl", "utils.jl"])
        @testset "$file" begin
            include(joinpath(@__DIR__, file))
        end
    end
end
