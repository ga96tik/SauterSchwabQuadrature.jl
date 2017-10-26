using CompScienceMeshes
export k2_cf

function k2_cf(u,v)
	û = [1-u[1], u[2]]
	v̂ = [1-v[1], v[2]]
	n1 = neighborhood(testchart, û)
	n2 = neighborhood(sourcechart, v̂)
	x = cartesian(n1)
	y = cartesian(n2)

	return(kernel(x,y)*jacobian(n1)*jacobian(n2))
end


#=function kernel2(x,y)
	#((testfunction(x))'*sourcefunction(y))*exp(-im*norm(x-y))/(4pi*norm(x-y))
	#can't store complex values in a matrix
	((testfunction(x))'*sourcefunction(y))*exp(norm(x-y))/(4pi*norm(x-y))
end=#

#=function kernel(x,y,pI,pII)

			return(((x-pI)'*(y-pII))*exp(-im*norm(x-y))/(4pi*norm(x-y)))
end=#
