using CompScienceMeshes
export kernel1
export kernel2

function kernel1(x,y)
	exp(-im*norm(x-y))/(4pi*norm(x-y))
	exp(norm(x-y))/(4pi*norm(x-y))
end

function kernel2(x,y)
	#((testfunction(x))'*sourcefunction(y))*exp(-im*norm(x-y))/(4pi*norm(x-y))
	#can't store complex values in a matrix
	((testfunction(x))'*sourcefunction(y))*exp(norm(x-y))/(4pi*norm(x-y))
end
