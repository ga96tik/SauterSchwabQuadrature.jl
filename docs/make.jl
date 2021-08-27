using Documenter, SauterSchwabQuadrature

makedocs(
    modules=[SauterSchwabQuadrature],
    doctest=false,
    sitename="SauterSchwabQuadrature",
)

# deploydocs(
#     deps = Deps.pip("pygments","mkdocs", "python-markdown-math"),
#     repo = "github.com/ga96tik/SauterSchwabQuadrature.jl.git",
#     julia  = "0.6",
#     osname = "linux")

deploydocs(
    # deps = Deps.pip("pygments","mkdocs", "python-markdown-math"),
    repo = "github.com/ga96tik/SauterSchwabQuadrature.jl.git",
    # julia  = "0.6",
    # osname = "linux",
)
