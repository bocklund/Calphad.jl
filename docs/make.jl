using Calphad
using Documenter

DocMeta.setdocmeta!(Calphad, :DocTestSetup, :(using Calphad); recursive=true)

makedocs(;
    modules=[Calphad],
    authors="Brandon Bocklund",
    repo="https://github.com/bocklund/Calphad.jl/blob/{commit}{path}#{line}",
    sitename="Calphad.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bocklund.github.io/Calphad.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bocklund/Calphad.jl",
    devbranch="main",
)
