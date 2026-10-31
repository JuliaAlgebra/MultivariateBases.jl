include("utils.jl")

for file in readdir(@__DIR__)
    if endswith(file, ".jl") &&
       !in(file, ["runtests.jl", "utils.jl"]) &&
       (isempty(ARGS) || file in ARGS)
        @testset "$file" begin
            include(joinpath(@__DIR__, file))
        end
    end
end
