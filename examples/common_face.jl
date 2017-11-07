using CompScienceMeshes
using SauterSchwabQuadrature



pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)

function integrand(x,y)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
	end

Sourcechart = Testchart = simplex(pI, pII, pIII)

cf = CommonFace(12)			#accuracy 12



result = sauterschwabintegral(Sourcechart, Testchart, integrand, cf)
println(result)
