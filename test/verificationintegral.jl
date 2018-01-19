include("parameterization.jl")

function verifintegral1(sourcechart, testchart, integrand, accuracy)

	global Kernel, p0, p1, p2

	Kernel = integrand

	p0 = Sourcechart.vertices[1]
	p1 = Sourcechart.vertices[2]
	p2 = Sourcechart.vertices[3]

	qps1 = quadpoints(Sourcechart, accuracy)

	path = simplex(point(0), point(1))
	qps2 = quadpoints(path, accuracy)

	result = sum(w*w2*w1*k(α,β,γ) for (β,w1) in qps2, (γ,w2) in qps2, (α,w) in qps1)

	return (result)

end



function verifintegral2(sourcechart, testchart, integrand, accuracy)

	qps1 = quadpoints(sourcechart, accuracy)
	qps2 = quadpoints(testchart, accuracy)

	result = sum(w2*w1*integrand(cartesian(x),cartesian(y))
				for (x,w2) in qps2, (y,w1) in qps1)

	return (result)

end
