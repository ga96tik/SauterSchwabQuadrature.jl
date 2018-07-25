using LinearAlgebra
using Test

using CompScienceMeshes
using SauterSchwabQuadrature
using StaticArrays
using BEAST

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pIV = point(5,1,-3)

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pII)

Accuracy = 12
ce = CommonEdge(Accuracy)

function integrand(x,y)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
end

function INTEGRAND(û,v̂)
	n1 = neighborhood(Testchart, û)
	n2 = neighborhood(Sourcechart, v̂)
	x = cartesian(n1)
	y = cartesian(n2)
	output = integrand(x,y)*jacobian(n1)*jacobian(n2)

return(output)
end

result = sauterschwab_parameterized(INTEGRAND, ce)-
		   verifintegral2(Sourcechart, Testchart, integrand, Accuracy)

@test norm(result) < 1.e-3



kernel(x,y) = 1/norm(cartesian(x)-cartesian(y))

t1 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.0, -0.980785, -0.19509],
    @SVector[0.0, -0.92388, -0.382683])
t2 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.373086, -0.881524, -0.289348],
    @SVector[0.294908, -0.944921, -0.141962])

@test indexin(t1.vertices, t2.vertices) == [1, nothing, nothing]

rt = BEAST.RTRefSpace{Float64}()
igd = generate_integrand_uv(kernel, rt, rt, t1, t2)

i5 = sauterschwab_parameterized(igd, CommonVertex(5))
i10 = sauterschwab_parameterized(igd, CommonVertex(10))

# brute numerical approach
q1 = quadpoints(t1, 10)
q2 = quadpoints(t2, 10)

M = N = numfunctions(rt)
iref = zero(i5)
for (x,w1) in q1
    f = rt(x)
    for (y,w2) in q2
        g = rt(y)
        G = kernel(x,y)
        ds = w1*w2
        global iref += SMatrix{M,N}([dot(f[i][1], G*g[j][1])*ds for i=1:M, j=1:N])
    end
end

@test i5  ≈ iref atol=1e-7
@test i10 ≈ iref atol=1e07
