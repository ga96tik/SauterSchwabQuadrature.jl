using CompScienceMeshes
using SauterSchwabQuadrature

include("parametrisation.jl")

function verifintegral1(Chart, Kernel)

	global kernel, chart, p0, p1, p2

	chart = Chart
	kernel = Kernel

	p0 = chart.vertices[1]
	p1 = chart.vertices[2]
	p2 = chart.vertices[3]

	acc = 12							#accuracy

	qps1 = quadpoints(chart, acc)

	path = simplex(point(0), point(1))
	qps2 = quadpoints(path, acc)

	result = sum(w*w2*w1*k(α,β,γ) for (β,w1) in qps2, (γ,w2) in qps2, (α,w) in qps1)

	return (result)

end



function verifintegral2(Sourcechart, Testchart, Kernel, ::Any)

	global kernel, testchart, sourcechart

	sourcechart = Sourcechart
	testchart = Testchart
	kernel = Kernel

	acc = 12								#accuracy

	qps1 = quadpoints(sourcechart, acc)
	qps2 = quadpoints(testchart, acc)

	result = sum(w2*w1*kernel(cartesian(x),cartesian(y))
				for (x,w2) in qps2, (y,w1) in qps1)

	return (result)

end
