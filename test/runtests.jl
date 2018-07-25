using SauterSchwabQuadrature

using LinearAlgebra
using Test

#include("test_ce_np_verification.jl")
#include("test_cf_np_verification.jl")
#include("test_cv_np_verification.jl")

include("test_cf_p_verification.jl")
include("test_ce_p_verification.jl")
include("test_cv_p_verification.jl")
include("test_pd_p_verification.jl")
