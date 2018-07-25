using LinearAlgebra
using CompScienceMeshes
using SauterSchwabQuadrature
using StaticArrays
using BEAST

kernel(x,y) = 1/norm(cartesian(x)-cartesian(y))

function integrand_fctr(kernel, testref, trialref, testel, trialel)

    function k3(u,v)
        out = @SMatrix zeros(3,3)

        x = neighborhood(testel,u)
        y = neighborhood(trialel,v)

        kernelval = kernel(x,y)
        f = testref(x)
        g = trialref(y)

        jx = jacobian(x)
        jy = jacobian(y)
        ds = jx*jy

        M = numfunctions(testref)
        N = numfunctions(trialref)

        return  SMatrix{M,N}([dot(f[i][1], kernelval*g[j][1])*ds for i=1:M, j=1:N])
    end

    return k3
end

fn = joinpath(dirname(@__FILE__), "sphere.in")
m = readmesh(fn)
X = raviartthomas(m)
x = refspace(X)

t1 = chart(m,cells(m)[1])
t2 = chart(m,cells(m)[100])

t1 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.0, -0.980785, -0.19509],
    @SVector[0.0, -0.92388, -0.382683])
t2 = simplex(
    @SVector[0.180878, -0.941848, -0.283207],
    @SVector[0.373086, -0.881524, -0.289348],
    @SVector[0.294908, -0.944921, -0.141962])

@assert indexin(t1.vertices, t2.vertices) == [1, nothing, nothing]

k = integrand_fctr(kernel, x, x, t1, t2)

i5 = sauterschwab_parameterized(k, CommonVertex(5))
i10 = sauterschwab_parameterized(k, CommonVertex(10))

# brute numerical approach
q1 = quadpoints(t1, 10)
q2 = quadpoints(t2, 10)

M = N = numfunctions(refspace(X))
iref = zero(i5)
for (x,w1) in q1
    f = refspace(X)(x)
    for (y,w2) in q2
        g = refspace(X)(y)
        G = kernel(x,y)
        ds = w1*w2
        iref += SMatrix{M,N}([dot(f[i][1], G*g[j][1])*ds for i=1:M, j=1:N])
    end
end
