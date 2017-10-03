using CompScienceMeshes

o, x, y, z = euclidianbasis(3)
p1 = point(0,0,0);
p2 = point(0,4,0);
p3 = point(5,1,0);
τ = simplex(p3,p2,p1)
volume(τ)

#c = neighborhood(τ,(10,2))
c = neighborhood(τ, (0.5,0.5))
parametric(c)
cartesian(c)
normal(c)
jacobian(c) #=== 2*volume(τ)=#


qps = quadpoints(τ, 6)
c = qps[1][1]
parametric(c)
cartesian(c)
normal(c)
jacobian(c) == 2*volume(τ)

f(p) = 1.0
sum(w*f(p) for (p,w) in qps)
volume(τ)

## Special case for integration over interval [0,1]
p1 = point(0.0)
p2 = point(1.0)

τ = simplex(p1,p2)
#volume(τ)

qps = quadpoints(τ,4)
#length(qps)

#f(p) = (x = cartesian(p)[1]; x^4)
function f(p)
  p = cartesian(p)[1]; p^2
end

sum(w*f(p) for (p,w) in qps)

## double integral over [0,1] × [0,1]
function f(p,q)
x = cartesian(p)[1]
y = cartesian(q)[1]
x*y
end


sum(wp*wq*f(p,q) for (p,wp) in qps, (q,wq) in qps)
