using CompScienceMeshes
using SauterSchwabQuadrature

include("parametrisation.jl")

function verifintegral(chart1, chart2, Kernel, ::common_faces)

	global kernel, sourcechart, testchart, p0, p1, p2

	sourcechart = chart1
	testchart = chart2
	kernel = Kernel

	p0 = sourcechart.vertices[1]
	p1 = sourcechart.vertices[2]
	p2 = sourcechart.vertices[3]


	qps1 = quadpoints(sourcechart, 12)

	path = simplex(point(0), point(1))
	qps2 = quadpoints(path, 12)

	result = sum(w*w2*w1*k(α,β,γ) for (β,w1) in qps2, (γ,w2) in qps2, (α,w) in qps1)

	return (result)

end
