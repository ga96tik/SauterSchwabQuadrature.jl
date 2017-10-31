using CompScienceMeshes

function k(α,β,γ)

			x̂ = cartesian(α)
			ξ = cartesian(β)[1]
			η = cartesian(γ)[1]

			ŷ = [ξ*(1-η), η]

			#jacobian of y--->ŷ = (1-y[2])

			chartI = simplex(p1,x̂,p0)
			chartII = simplex(p2,x̂,p1)
			chartIII = simplex(p0,x̂,p2)

			n1 = neighborhood(chartI, ŷ)
			n2 = neighborhood(chartII, ŷ)
			n3 = neighborhood(chartIII, ŷ)


			yI = cartesian(n1)
			yII = cartesian(n2)
			yIII = cartesian(n3)

			return(Kernel(x̂,yI)*jacobian(n1)*(1-η) +
					Kernel(x̂,yII)*jacobian(n2)*(1-η) +
					Kernel(x̂,yIII)*jacobian(n3)*(1-η))
end
