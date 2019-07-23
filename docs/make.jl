using Documenter, MultivariateBases

makedocs(
    sitename = "MultivariateBases",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Introduction" => "index.md",
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [MultivariateBases]
)

deploydocs(
    repo   = "github.com/JuliaAlgebra/MultivariateBases.jl.git",
)
