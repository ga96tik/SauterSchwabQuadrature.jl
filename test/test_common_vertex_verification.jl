using CompScienceMeshes
using Base.Test
using SauterSchwabQuadrature

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pIV = point(5,1,-3)
pV = point(0,0,0)

function integrand(x,y)
			pI = point(1,5,3)
			pII = point(2,5,3)
			return(((x-pI)'*(y-pV))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
end

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pV)

Accuracy = 12
cv = CommonVertex(Accuracy)

result = sauterschwab_nonparameterized(Sourcechart, Testchart, integrand, cv) -
		verifintegral2(Sourcechart, Testchart, integrand, Accuracy)

@test norm(result) < 1.e-3
