using CompScienceMeshes
using Base.Test
using SauterSchwabQuadrature

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
Sourcechart = Testchart = simplex(pI, pII, pIII)

function integrand(x,y)
			pI = point(1,5,3)
			pII = point(2,5,3)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end

cf = CommonFace(12)

result = verifintegral1(Sourcechart, Testchart, integrand, cf) -
	sauterschwabintegral(Sourcechart, Testchart, integrand, cf)

@test norm(result) < 1.e-3
