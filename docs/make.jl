using AMGCLWrap
using Documenter

DocMeta.setdocmeta!(AMGCLWrap, :DocTestSetup, :(using AMGCLWrap); recursive=true)

makedocs(;
    modules=[AMGCLWrap],
    authors="Jürgen Fuhrmann <juergen-fuhrmann@web.de> and contributors",
    repo="https://github.com/j-fu/AMGCLWrap.jl/blob/{commit}{path}#{line}",
    sitename="AMGCLWrap.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://j-fu.github.io/AMGCLWrap.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/j-fu/AMGCLWrap.jl",
    devbranch="main",
)
