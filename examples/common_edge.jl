using CompScienceMeshes
using SauterSchwabQuadrature



pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pIV = point(5,1,-3)

function integrand(x,y)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
end

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pII)

ce = CommonEdge(12)			#accuracy 12



result = sauterschwabintegral(Sourcechart, Testchart, integrand, ce)
println(result)
