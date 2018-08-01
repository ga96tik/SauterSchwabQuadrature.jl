using LinearAlgebra
using CompScienceMeshes
using BEAST
using SauterSchwabQuadrature
using StaticArrays

t1 = simplex(
	@SVector[0.180878, -0.941848, -0.283207],
 	@SVector[0.0, -0.980785, -0.19509],
	@SVector[0.0, -0.92388, -0.382683],)
t2 = simplex(
	@SVector[0.180878, -0.941848, -0.283207],
	@SVector[0.0, -0.92388, -0.382683],
	@SVector[0.158174, -0.881178, -0.44554],)

rt = BEAST.RTRefSpace{Float64}()
kernel = (x,y) -> 1/norm(cartesian(x)-cartesian(y))
igd = generate_integrand_uv(kernel, rt, rt, t1, t2)

u = (1/3,1/3,)
v = (1/3,1/3,)

igd(u,v)

i15 = sauterschwab_parameterized(igd, CommonEdge(15))
