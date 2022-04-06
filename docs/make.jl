using AdaptiveHierarchicalRegularBinning
using Documenter

DocMeta.setdocmeta!(AdaptiveHierarchicalRegularBinning, :DocTestSetup, :(using AdaptiveHierarchicalRegularBinning); recursive=true)

makedocs(;
    modules=[AdaptiveHierarchicalRegularBinning],
    authors="Nikolaos Pitsianis <nikos.pitsianis@eng.auth.gr> and contributors",
    repo="https://github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl/blob/{commit}{path}#{line}",
    sitename="AdaptiveHierarchicalRegularBinning.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pitsianis.github.io/AdaptiveHierarchicalRegularBinning.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pitsianis/AdaptiveHierarchicalRegularBinning.jl",
    devbranch="main",
)
