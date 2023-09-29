using Documenter, SauterSchwabQuadrature

DocMeta.setdocmeta!(SauterSchwabQuadrature, :DocTestSetup, :(using SauterSchwabQuadrature); recursive=true)

makedocs(;
    modules=[SauterSchwabQuadrature],
    sitename="SauterSchwabQuadrature.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://ga96tik.github.io/SauterSchwabQuadrature.jl", edit_link="master", assets=String[]
    ),
    pages=[
        "Introduction" => "index.md",
        "Details" => "details.md",
        "Manual" => "manual.md",
        "API Reference" => "apiref.md",
    ],
)

deploydocs(; repo="github.com/ga96tik/SauterSchwabQuadrature.jl.git", devbranch="master")
