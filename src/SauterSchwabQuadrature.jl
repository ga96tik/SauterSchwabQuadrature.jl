module SauterSchwabQuadrature

using LinearAlgebra

using CompScienceMeshes
using StaticArrays

export sauterschwabintegral, sauterschwab_nonparameterized, sauterschwab_parameterized
export CommonFace, CommonEdge, CommonVertex, PositiveDistance
export generate_integrand_uv

include("k3.jl")
#include("parameterization.jl")
include("sauterschwabintegral.jl")

include("parametric_kernel_generator.jl")

end
