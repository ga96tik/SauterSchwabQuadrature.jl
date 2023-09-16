module SauterSchwabQuadrature


# -------- used packages
using LinearAlgebra
using StaticArrays



# -------- exportet parts
# types
export SauterSchwabStrategy
export CommonFace, CommonEdge, CommonVertex, PositiveDistance
export CommonFaceQuad, CommonEdgeQuad, CommonVertexQuad

# functions
export sauterschwab_parameterized



# -------- included files
include("sauterschwabintegral.jl")
include("pulled_back_integrands.jl")
include("reorder_vertices.jl")


end
