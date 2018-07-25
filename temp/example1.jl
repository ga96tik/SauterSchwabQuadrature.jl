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

        return  SMatrix([dot(f[i][1], kernelval*g[j][1])*ds for i=1:M, j=1:N])
    end

    return k3
end

fn = joinpath(dirname(@__FILE__), "sphere.in")
m = readmesh(fn)
X = raviartthomas(m)
x = refspace(X)

t1 = chart(m,cells(m)[1])
t2 = chart(m,cells(m)[100])

k = integrand_fctr(kernel, x, x, t1, t2)
