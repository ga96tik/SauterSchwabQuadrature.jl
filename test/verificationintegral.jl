using CompScienceMeshes
using SauterSchwabQuadrature

include("parametrisation.jl")

function verifintegral1(sourcechart, testchart, kernel, ::Any)

	global Kernel, Sourcechart, Testchart, p0, p1, p2

	Sourcechart = sourcechart
	Testchart = testchart
	Kernel = kernel

	p0 = Sourcechart.vertices[1]
	p1 = Sourcechart.vertices[2]
	p2 = Sourcechart.vertices[3]

	acc = 12							#accuracy

	qps1 = quadpoints(Sourcechart, acc)

	path = simplex(point(0), point(1))
	qps2 = quadpoints(path, acc)

	result = sum(w*w2*w1*k(α,β,γ) for (β,w1) in qps2, (γ,w2) in qps2, (α,w) in qps1)

	return (result)

end



function verifintegral2(sourcechart, testchart, kernel, accuracy::Any)

	global Kernel, Testchart, Sourcechart

	Sourcechart = sourcechart
	Testchart = testchart
	Kernel = kernel

	Accuracy = accuracy

	qps1 = quadpoints(Sourcechart, Accuracy.acc)
	qps2 = quadpoints(Testchart, Accuracy.acc)

	result = sum(w2*w1*Kernel(cartesian(x),cartesian(y))
				for (x,w2) in qps2, (y,w1) in qps1)

	return (result)

end
