export sauterschwabintegral, sauterschwab_nonparameterized, sauterschwab_parameterized
export CommonFace, CommonEdge, CommonVertex, PositiveDistance

include(Pkg.dir("SauterSchwabQuadrature","test","verificationintegral.jl"))


struct CommonFace <: Any		acc::Int64		end
struct CommonEdge <:Any			acc::Int64		end
struct CommonVertex <:Any		acc::Int64		end
struct PositiveDistance <:Any	acc::Int64		end





function sauterschwabintegral(sourcechart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
	 							testchart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
								 integrand, accuracy::Int64, accuracy_pd::Int64)

	index_equal_src = findin(sourcechart.vertices, testchart.vertices)
	index_equal_tst = findin(testchart.vertices, sourcechart.vertices)

	m = length(index_equal_src)


	if m == 0
		Sourcechart = sourcechart
		Testchart = testchart
		method = PositiveDistance(accuracy_pd)

	elseif m == 1
		SOURCECHART = Array{Any}(3)
		TESTCHART = Array{Any}(3)

		SOURCECHART[2] = sourcechart.vertices[2]
		TESTCHART[2] = testchart.vertices[2]
		SOURCECHART[3] = sourcechart.vertices[3]
		TESTCHART[3] = testchart.vertices[3]

		SOURCECHART[1] = sourcechart.vertices[index_equal_src[1]]
		TESTCHART[1] = sourcechart.vertices[index_equal_src[1]]
		SOURCECHART[index_equal_src[1]] = sourcechart.vertices[1]
		TESTCHART[index_equal_tst[1]] = testchart.vertices[1]
		Sourcechart = simplex(SOURCECHART[1], SOURCECHART[2], SOURCECHART[3])
		Testchart = simplex(TESTCHART[1], TESTCHART[2], TESTCHART[3])
		method = CommonVertex(accuracy)

	elseif m == 2
		SOURCECHART = Array{Any}(3)
		TESTCHART = Array{Any}(3)

		SOURCECHART[2] = sourcechart.vertices[2]
		TESTCHART[2] = testchart.vertices[2]

		SOURCECHART[1] = TESTCHART[1] = sourcechart.vertices[index_equal_src[1]]
		SOURCECHART[3] = TESTCHART[3] = sourcechart.vertices[index_equal_src[2]]
		SOURCECHART[index_equal_src[1]] = sourcechart.vertices[1]
		SOURCECHART[index_equal_src[2]] = sourcechart.vertices[3]
		TESTCHART[index_equal_tst[1]] = testchart.vertices[1]
		TESTCHART[index_equal_tst[2]] = testchart.vertices[3]
		Sourcechart = simplex(SOURCECHART[1], SOURCECHART[2], SOURCECHART[3])
		Testchart = simplex(TESTCHART[1], TESTCHART[2], TESTCHART[3])
		method = CommonEdge(accuracy)

	else
		Sourcechart = Testchart = sourcechart
		method = CommonFace(accuracy)
	end


	sauterschwab_nonparameterized(Sourcechart, Testchart, integrand, method)

end






function sauterschwab_nonparameterized(sourcechart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
	 									testchart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
										integrand, method::CommonFace)

	global Sourcechart, Testchart, Kernel
	Sourcechart = sourcechart			#sourcetriangle t
	Testchart = testchart				#testtriangle τ
	Kernel = integrand

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, method.acc)

		result = sum(w1*w2*w3*w4*k3_cf(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwab_nonparameterized(sourcechart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
	 									testchart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
										integrand, method::CommonEdge)

#Sourcechart and Testchart must have been generated in the following manner:
#Sourcechart = simplex(a,b,c);	Testchart = simplex(a,b',c)
#First and third entry of both charts must be equal

	global Kernel, Sourcechart, Testchart
	Kernel = integrand
	Sourcechart = sourcechart		#sourcetriangle t
	Testchart = testchart			#testtriangle τ

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, method.acc)

	result = sum(w1*w2*w3*w4*k3_ce(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwab_nonparameterized(sourcechart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
	 									testchart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
										integrand, method::CommonVertex)

#Sourcechart and Testchart must have been generated in the following manner:
#Sourcechart = simplex(a,b,c);	Testchart = simplex(a,b',c')
#First entry of both charts must be equal

	global Kernel, Sourcechart, Testchart
	Kernel = integrand
	Sourcechart = sourcechart		#sourcetriangle t
	Testchart = testchart			#testtriangle τ

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, method.acc)

	result = sum(w1*w2*w3*w4*k3_cv(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwab_nonparameterized(sourcechart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
	 									testchart::CompScienceMeshes.Simplex{3,2,1,3,Float64},
										integrand, method::PositiveDistance)

	result = verifintegral2(sourcechart, testchart, integrand, method.acc)

	return(result)
end

















function sauterschwab_parameterized(sourcechart, testchart, integrand, method::CommonFace)

	global  Kernel
	Kernel = integrand

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, method.acc)

		result = sum(w1*w2*w3*w4*k3p_cf(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwab_parameterized(sourcechart, testchart, integrand, method::CommonEdge)

	global Kernel
	Kernel = integrand

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, method.acc)

	result = sum(w1*w2*w3*w4*k3p_ce(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwab_parameterized(sourcechart, testchart, integrand, method::CommonVertex)

	global Kernel
	Kernel = integrand

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, method.acc)

	result = sum(w1*w2*w3*w4*k3p_cv(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end





function sauterschwab_parameterized(sourcechart, testchart, integrand, method::PositiveDistance)

	global Kernel
	Kernel = integrand

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, method.acc)

	result = sum(w1*w2*w3*w4*k3p_pd(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)
end
