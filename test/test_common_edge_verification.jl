using CompScienceMeshes
using Base.Test
using SauterSchwabQuadrature

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pIV = point(5,1,-3)

function kernel(x,y)
			pI = point(1,5,3)
			pII = point(2,5,3)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pII)
ce = CommonEdge(12)

result = sauterschwabintegral(Sourcechart, Testchart, kernel, ce) -
		verifintegral2(Sourcechart, Testchart, kernel, ce)

@test norm(result) < 1.e-3
