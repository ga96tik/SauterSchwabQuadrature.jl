using CompScienceMeshes
using SauterSchwabQuadrature
using LinearAlgebra 

ô = point(0,0,0)
x̂ = point(1,0,0)
ŷ = point(0,1,0)
ẑ = point(0,0,1)

triangle1 = simplex(2x̂, 3ŷ, ô)
triangle2 = simplex(2x̂, 4x̂+3ŷ+ẑ, ô)


function fxy(p,q)
    x = cartesian(p)
    y = cartesian(q)
    n = normal(p)
    R = norm(x-y)
    1/R
    # dot(n, x-y)/R^3
end

function fuv(u,v)
    p = neighborhood(triangle1, u)
    q = neighborhood(triangle2, v)
    J = jacobian(p) * jacobian(q)
    return fxy(p, q) * J
end

nmax = 30
results_ce = zeros(nmax)
results_pd = zeros(nmax)
errors = zeros(nmax)
errors[1] = 1
for n in 1:nmax
    ce = CommonEdge(SauterSchwabQuadrature._legendre(n,0.0,1.0))
    pd = PositiveDistance(SauterSchwabQuadrature._legendre(n,0.0,1.0))
    results_ce[n] = sauterschwab_parameterized(fuv, ce)
    results_pd[n] = sauterschwab_parameterized(fuv, pd)
    n > 1 && (errors[n] = abs(results_ce[n] - results_ce[n-1]))
    @show errors[n]
end

exact = results_ce[end]
@assert abs(errors[end]) < 1e-5

using Plots
plot(title="convergence SS quadrature")
plot!(log10.(abs.(exact .- results_ce)), label="common edge [A,P1,B]-[A,P2,B]")
plot!(log10.(abs.(exact .- results_pd)), label="positive distance")
# plot!(log10.(abs.(exact .- results3)))
# plot!(log10.(abs.(exact .+ resultspd)))