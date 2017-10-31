using CompScienceMeshes
export sauterschwabintegral, CommonFace, CommonEdge, CommonVertex, PositiveDistance
include(Pkg.dir("SauterSchwabQuadrature","test","verificationintegral.jl"))

struct CommonFace <: Any		acc::Int64		end
struct CommonEdge <: Any 		acc::Int64		end
struct CommonVertex <:Any		acc::Int64		end
struct PositiveDistance <:Any	acc::Int64		end






function sauterschwabintegral(sourcechart, testchart, kernel, accuracy::CommonFace)

	global Sourcechart, Testchart, Kernel
	Sourcechart = sourcechart
	Testchart = testchart
	Kernel = kernel

	if Testchart == Sourcechart
		nothing
	else
		println("Error: Testchart \u2260 Sourcechart")
		quit()
	end

	Accuracy = accuracy
	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, Accuracy.acc)

		result = sum(w1*w2*w3*w4*k3_cf(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwabintegral(sourcechart, testchart, kernel, accuracy::CommonEdge)

#Sourcechart and Testchart must have been generated in the following manner:
#Sourcechart = simplex(a,b,c);	Testchart = simplex(a,b',c)
#First and third entry of both charts must be equal

	global Kernel, Sourcechart, Testchart
	Kernel = kernel
	Sourcechart = sourcechart		#sourcetriangle t
	Testchart = testchart			#testtriangle τ

	Accuracy = accuracy
	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, Accuracy.acc)

	result = sum(w1*w2*w3*w4*k3_ce(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwabintegral(sourcechart, testchart, kernel, accuracy::CommonVertex)

#Sourcechart and Testchart must have been generated in the following manner:
#Sourcechart = simplex(a,b,c);	Testchart = simplex(a,b',c')
#First entry of both charts must be equal

	global Kernel, Sourcechart, Testchart
	Kernel = kernel
	Sourcechart = sourcechart		#sourcetriangle t
	Testchart = testchart			#testtriangle τ

	Accuracy = accuracy
	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, Accuracy.acc)

	result = sum(w1*w2*w3*w4*k3_cv(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwabintegral(sourcechart, testchart, kernel, accuracy::PositiveDistance)

#The intersection of testtriangle and sourcetriangle must be 0!

	global Kernel, Sourcechart, Testchart
	Kernel = kernel
	Sourcechart = sourcechart		#sourcetriangle t
	Testchart = testchart			#testtriangle τ

	Accuracy = accuracy

	result = verifintegral2(Sourcechart, Testchart, Kernel, Accuracy)

	return(result)
end
