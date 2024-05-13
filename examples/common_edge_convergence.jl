using CompScienceMeshes
using BEAST
using SauterSchwabQuadrature
using LinearAlgebra

τ = simplex(
    point(-0.9423169199664047, 0.32954812598003336, -0.011695233370427325),
    point(-0.9551397067786591, 0.2847288579921698, -0.00779682224695155),
    point(-0.9238795325109631, 0.382683432365871, 3.0782976471689056e-18),
)

σ = simplex(
    point(-0.9807852804031678, 0.19509032201644239, -7.632783294297951e-17),
    point(-0.9238795325109631, 0.38268343236587105, -1.7210136298535816e-16),
    point(-0.9551397067786591, 0.2847288579921698, -0.00779682224695155),
)

τ2 = simplex(
    point(-0.9238795325109631, 0.382683432365871, 3.0782976471689056e-18),
    point(-0.9423169199664047, 0.32954812598003336, -0.011695233370427325),
    point(-0.9551397067786591, 0.2847288579921698, -0.00779682224695155),
)

σ2 = simplex(
    point(-0.9238795325109631, 0.38268343236587105, -1.7210136298535816e-16),
    point(-0.9807852804031678, 0.19509032201644239, -7.632783294297951e-17),
    point(-0.9551397067786591, 0.2847288579921698, -0.00779682224695155),
)

τ3 = simplex(
    point(-0.9551397067786591, 0.2847288579921698, -0.00779682224695155),
    point(-0.9423169199664047, 0.32954812598003336, -0.011695233370427325),
    point(-0.9238795325109631, 0.382683432365871, 3.0782976471689056e-18),
)

σ3 = simplex(
    point(-0.9551397067786591, 0.2847288579921698, -0.00779682224695155),
    point(-0.9807852804031678, 0.19509032201644239, -7.632783294297951e-17),
    point(-0.9238795325109631, 0.38268343236587105, -1.7210136298535816e-16),
)


ϕ = BEAST.RTRefSpace{Float64}()

function integrand(x, y)
    R = norm(x - y)
    κ = 1.0
    exp(-im * κ * R) / R / 4 / pi - 1 / R / 4 / pi
end

function INTEGRAND(û, v̂)
    n1 = neighborhood(τ2, û)
    n2 = neighborhood(σ2, v̂)
    x = cartesian(n1)
    y = cartesian(n2)
    fx = ϕ(n1)[1]
    gy = ϕ(n2)[2]
    return integrand(x, y) * jacobian(n1) * jacobian(n2) * (dot(fx.value, gy.value) - fx.divergence * gy.divergence)
end

results = []
for nodes in 1:20
    ce = SauterSchwabQuadrature.CommonEdge(SauterSchwabQuadrature._legendre(nodes, 0.0, 1.0))
    result = SauterSchwabQuadrature.sauterschwab_parameterized(INTEGRAND, ce)
    push!(results, result)
end

errors = abs.((results .- results[end]) / results[end])
using Plots: Plots
Plots.plot(log10.(abs.(errors[1:(end - 1)])))


p1 = point(0.0, 0.0, 0.0) # the same for both quads
p2 = point(2.0, 0.0, 0.0) # the same for both quads
p3 = point(2.0, 2.0, 0.0)
p4 = point(0.0, 2.0, 0.0)

p3 = point(2.0, -2.0, 0.0) # choose points of second quad such that parametrizations align as required in [1]
p4 = point(0.0, -2.0, 0.0)

using StaticArrays
struct Quadrilateral
    p1::SVector{3}
    p2::SVector{3}
    p3::SVector{3}
    p4::SVector{3}
end

function jacobianDet(quad::Quadrilateral, u)

    aux = quad.p3 - quad.p4 + quad.p1 - quad.p2

    ∂ru = quad.p2 - quad.p1 + u[2] * aux
    ∂rv = quad.p4 - quad.p1 + u[1] * aux

    D = (∂ru[2] * ∂rv[3] - ∂ru[3] * ∂rv[2])^2 + (∂ru[3] * ∂rv[1] - ∂ru[1] * ∂rv[3])^2 + (∂ru[1] * ∂rv[2] - ∂ru[2] * ∂rv[1])^2

    return sqrt(D)
end

# parametrize planar quadrilateral with u, v ∈ [0,1]
function (quad::Quadrilateral)(u, v)
    return quad.p1 + u * (quad.p2 - quad.p1) + v * (quad.p4 - quad.p1) + u * v * (quad.p3 - quad.p4 + quad.p1 - quad.p2) # see, e.g., [1] page 187
end

struct singularKernel
    quad1::Quadrilateral
    quad2::Quadrilateral
end


function (sK::singularKernel)(u, v)

    x = sK.quad1(u...)
    y = sK.quad2(v...)

    return jacobianDet(sK.quad1, u) * jacobianDet(sK.quad2, v) / norm(x - y)
end

q1 = Quadrilateral(p1, p2, p3, p4)
q2 = Quadrilateral(p1, p2, p3, p4)

sK = singularKernel(q1, q2)

results = []
for nodes in 1:20
    ce = SauterSchwabQuadrature.CommonEdge(SauterSchwabQuadrature._legendre(nodes, 0.0, 1.0))
    result = SauterSchwabQuadrature.sauterschwab_parameterized(sK, ce)
    push!(results, result)
end
errors = abs.((results .- results[end]) / results[end])
Plots.plot(log10.(abs.(errors[1:(end - 1)])))
