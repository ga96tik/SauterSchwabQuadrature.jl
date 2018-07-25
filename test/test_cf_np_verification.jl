using CompScienceMeshes
using Test
using SauterSchwabQuadrature

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
Sourcechart = Testchart = simplex(pI, pII, pIII)

function integrand(x,y)
	return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
end

Accuracy = 12
cf = CommonFace(Accuracy)

result = verifintegral1(Sourcechart, Testchart, integrand, Accuracy) -
	sauterschwab_nonparameterized(Sourcechart, Testchart, integrand, cf)

@test norm(result) < 1.e-3
