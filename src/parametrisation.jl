using CompScienceMeshes
export k2

function k2(u,v)
	û = [1-u[1], u[2]]
	v̂ = [1-v[1], v[2]]
	n1 = neighborhood(testchart, û)
	n2 = neighborhood(sourcechart, v̂)
	x = cartesian(n1)
	y = cartesian(n2)

	return(kernel(x,y)*jacobian(n1)*jacobian(n2))
end
