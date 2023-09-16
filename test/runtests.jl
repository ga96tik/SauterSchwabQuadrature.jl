using SauterSchwabQuadrature
using Test

using LinearAlgebra
using StaticArrays
using CompScienceMeshes

# --- testsets
@testset "Testing SauterSchwabQuadrature functionality" begin

    include("parametric_kernel_generator.jl")
    include("local_space.jl")
    include("numquad.jl")
    include("verificationintegral.jl")

    @testset "Triangular " begin
        include("test_reordering.jl")
        include("test_cf_tr.jl")
        include("test_ce_tr.jl")
        include("test_cv_tr.jl")
        include("test_pd_tr.jl")
    end

    @testset "Quadrilateral " begin
        include("test_cf_quad.jl")
    end

end