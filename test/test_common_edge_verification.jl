using CompScienceMeshes
using Base.Test
using SauterSchwabQuadrature

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pIV = point(5,1,-3)

function integrand(x,y)
			pI = point(1,5,3)
			pII = point(2,5,3)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
end

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pII)

Accuracy = 12
ce = CommonEdge(Accuracy)

result = verifintegral2(Sourcechart, Testchart, integrand, Accuracy) -
			sauterschwab_nonparameterized(Sourcechart, Testchart, integrand, ce)

@test norm(result) < 1.e-3
