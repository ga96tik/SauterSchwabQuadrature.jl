using CompScienceMeshes
export sauterschwabintegral, common_edge, common_vertex, positive_distance
include(Pkg.dir("SauterSchwabQuadrature","test","verificationintegral.jl"))

type common_edge <: Any 		end
type common_vertex <:Any		end
type positive_distance <:Any	end






function sauterschwabintegral(Chart, Kernel)

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, 12)          #accuracy 12

	global sourcechart, testchart, kernel
	sourcechart = testchart = Chart
	kernel = Kernel

	result = sum(w1*w2*w3*w4*k3_cf(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwabintegral(Sourcechart, Testchart, Kernel, ::common_edge)

#Sourcechart and Testchart must have been generated in the following manner:
#Sourcechart = simplex(a,b,c);	Testchart = simplex(a,b',c)
#First and third entry of both charts must be equal

	global kernel, sourcechart, testchart
	kernel = Kernel
	sourcechart = Sourcechart		#sourcetriangle t
	testchart = Testchart			#testtriangle τ

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, 12)          #accuracy 12

	result = sum(w1*w2*w3*w4*k3_ce(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwabintegral(Sourcechart, Testchart, Kernel, ::common_vertex)

#Sourcechart and Testchart must have been generated in the following manner:
#Sourcechart = simplex(a,b,c);	Testchart = simplex(a,b',c')
#First entry of both charts must be equal

	global kernel, sourcechart, testchart
	kernel = Kernel
	sourcechart = Sourcechart		#sourcetriangle t
	testchart = Testchart			#testtriangle τ

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, 12)          #accuracy 12

	result = sum(w1*w2*w3*w4*k3_cv(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwabintegral(Sourcechart, Testchart, Kernel, ::positive_distance)

#The intersection of testtriangle and sourcetriangle must be 0!

	global kernel, sourcechart, testchart
	kernel = Kernel
	sourcechart = Sourcechart		#sourcetriangle t
	testchart = Testchart			#testtriangle τ

	result = verifintegral2(sourcechart, testchart, kernel, nothing)

	return(result)
end
