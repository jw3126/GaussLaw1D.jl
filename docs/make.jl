using Documenter, GaussLaw1D

makedocs(;
    modules=[GaussLaw1D],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jw3126/GaussLaw1D.jl/blob/{commit}{path}#L{line}",
    sitename="GaussLaw1D.jl",
    authors="Jan Weidner",
    assets=String[],
)
