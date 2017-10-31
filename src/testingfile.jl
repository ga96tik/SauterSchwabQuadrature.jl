using CompScienceMeshes
using SauterSchwabQuadrature

global pI, pII, pIII, pIV, pV, pVI, pVII, pVIII

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pIV = point(5,1,-3)
pV = point(0,0,0)
pVI = point(10,11,12)
pVII = point(10,11,13)
pVIII = point(11,11,12)




#common_face
function kernel(x,y)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end

Sourcechart = Testchart = simplex(pI, pII, pIII)

cf = CommonFace(12)

result = sauterschwabintegral(Sourcechart, Testchart, kernel, cf)






#common_edge
function kernel(x,y)
			return(((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pII)

ce = CommonEdge(12)

result = sauterschwabintegral(Sourcechart, Testchart, kernel, ce)







#common_vertex
function kernel(x,y)
			return(((x-pI)'*(y-pV))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pV)

cv = CommonVertex(12)

result = sauterschwabintegral(Sourcechart, Testchart, kernel, cv)






#positive_distance
function kernel(x,y)
			return(((x-pI)'*(y-pVII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))	#k=1
end

#Sourcechart and Testchart are chosen in a way, so that the intersection of both charts is 0.
Sourcechart = simplex(pI,pII,pIII)
Testchart = simplex(pVI,pVII,pVIII)

pd = PositiveDistance(12)

result = sauterschwabintegral(Sourcechart, Testchart, kernel, pd)
