using CompScienceMeshes
export sauterschwabintegral!

function sauterschwabintegral!(sourcecellid, testcellid, mesh)
	srccellvert = mesh.geo.faces[sourcecellid]
	tstcellvert = mesh.geo.faces[testcellid]
	zlocal = zeros(3,3)



	#figuring out of integration strategy
	strategy = 0
		for i = 1:3
			for j = 1:3
				if srccellvert[i] == tstcellvert[j]
					strategy = strategy+1
				else
				end
			end
		end



	if strategy == 0
		#positive distance
		println("0")

	elseif strategy == 1
		#common vertex
		println("1")

	elseif strategy == 2
		#common edge
		println("2")

	else
		for i = 1:3
			for j = 1:3

				global pI, pII, p0, p1, p2

				pI = mesh.geo.vertices[mesh.geo.faces[sourcecellid][j]]
				pII = mesh.geo.vertices[mesh.geo.faces[testcellid][i]]

				sourcecellid == testcellid
				id = sourcecellid
				vertice_ids = mesh.geo.faces[id]
				p0 = mesh.geo.vertices[vertice_ids[1]]
				p1 = mesh.geo.vertices[vertice_ids[2]]
				p2 = mesh.geo.vertices[vertice_ids[3]]

				start = point(0)
				stop = point(1)
				path = simplex(start, stop)
				qps = quadpoints(path, 2)          #accuracy 2

				zlocal[i,j] = sum(w1*w2*w3*w4*k3_2(η1, η2, η3, ξ)
					for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps)

			end
		end
	end

	println(zlocal)

end

#=function sauterschwabintegral!(sourcecellid, testcellid, mesh)
id = sourcecellid
vertice_ids = mesh.geo.faces[id]

global p0, p1, p2

p0 = mesh.geo.vertices[vertice_ids[1]]
p1 = mesh.geo.vertices[vertice_ids[2]]
p2 = mesh.geo.vertices[vertice_ids[3]]

start = point(0)
stop = point(1)
path = simplex(start, stop)
qps = quadpoints(path, 2)          #accuracy 2

return(sum(w1*w2*w3*w4*k31(η1, η2, η3, ξ) for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps))

end


function sauterschwabintegral!(sourcecellid, testcellid, mesh,
sourcefunctionid, testfunctionid)

global pI, pII, p0, p1, p2

pI = mesh.geo.vertices[mesh.geo.faces[sourcecellid][sourcefunctionid]]
pII = mesh.geo.vertices[mesh.geo.faces[testcellid][testfunctionid]]

sourcecellid == testcellid
id = sourcecellid
vertice_ids = mesh.geo.faces[id]
p0 = mesh.geo.vertices[vertice_ids[1]]
p1 = mesh.geo.vertices[vertice_ids[2]]
p2 = mesh.geo.vertices[vertice_ids[3]]

start = point(0)
stop = point(1)
path = simplex(start, stop)
qps = quadpoints(path, 2)          #accuracy 2

return(sum(w1*w2*w3*w4*k32(η1, η2, η3, ξ) for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ, w4) in qps))

end=#
