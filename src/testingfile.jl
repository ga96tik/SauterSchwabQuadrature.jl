using CompScienceMeshes
using BEAST
using SauterSchwabQuadrature

mesh = readmesh(Pkg.dir("SauterSchwabQuadrature", "src", "testingsphere.in"))
mesh = raviartthomas(mesh)

sauterschwabintegral!(1, 1, mesh)
sauterschwabintegral!(2, 2, mesh)
