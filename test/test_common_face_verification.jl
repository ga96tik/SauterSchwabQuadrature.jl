using CompScienceMeshes
using Base.Test
using SauterSchwabQuadrature

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
chart = simplex(pI, pII, pIII)

function kernel(x,y)
			pI = point(1,5,3)
			pII = point(2,5,3)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end

result = verifintegral1(chart, kernel) -
	sauterschwabintegral(chart, kernel)

@test norm(result) < 1.e-3
