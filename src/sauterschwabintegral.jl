using CompScienceMeshes
export sauterschwabintegral, common_faces

type common_faces <: Any
end

function sauterschwabintegral(chart1, chart2, Kernel, ::common_faces)

	start = point(0)
	stop = point(1)
	path = simplex(start, stop)
	qps = quadpoints(path, 12)          #accuracy 12

	global sourcechart, testchart, kernel
	sourcechart = chart1
	testchart = chart2
	kernel = Kernel


	result = sum(w1*w2*w3*w4*k3_cf(η1, η2, η3, ξ)
	for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

	return(result)

end
