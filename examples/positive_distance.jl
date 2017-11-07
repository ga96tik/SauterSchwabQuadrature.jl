using CompScienceMeshes
using SauterSchwabQuadrature



pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pVI = point(10,11,12)
pVII = point(10,11,13)
pVIII = point(11,11,12)

function integrand(x,y)
			return(((x-pI)'*(y-pVII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
end

Sourcechart = simplex(pI,pII,pIII)
Testchart = simplex(pVI,pVII,pVIII)

pd = PositiveDistance(12)			#accuracy 12



result = sauterschwabintegral(Sourcechart, Testchart, integrand, pd)
println(result)
