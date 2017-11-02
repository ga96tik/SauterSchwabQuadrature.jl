# #push!(LOAD_PATH,"C::/Users/Engelbert/.julia/SauterSchwabQuadrature/docs/src/")
# using Documenter, SauterSchwabQuadrature
#
# makedocs(modules=[SauterSchwabQuadrature],
#         doctest=false)
#
# deploydocs(deps = Deps.pip("pygments","mkdocs", "python-markdown-math"),
#     repo = "github.com/ga96tik/SauterSchwabQuadrature.git",
#     julia  = "0.6",
#     osname = "linux")


using Documenter, SauterSchwabQuadrature

makedocs(clean=false)
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/ga96tik/SauterSchwabQuadrature.jl.git",
    julia  = "0.6"
)
