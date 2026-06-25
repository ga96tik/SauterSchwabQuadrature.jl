using FastGaussQuadrature

abstract type SauterSchwabStrategy end

struct CommonFace{A} <: SauterSchwabStrategy
    qps::A
end
struct CommonEdge{A} <: SauterSchwabStrategy
    qps::A
end
struct CommonVertex{A} <: SauterSchwabStrategy
    qps::A
end
struct PositiveDistance{A} <: SauterSchwabStrategy
    qps::A
end


struct CommonFaceQuad{A} <: SauterSchwabStrategy
    qps::A
end
struct CommonEdgeQuad{A} <: SauterSchwabStrategy
    qps::A
end
struct CommonVertexQuad{A} <: SauterSchwabStrategy
    qps::A
end


function _legendre(n, a, b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b - a) / 2
    x = (x .+ 1) / 2 * (b - a) .+ a
    collect(zip(x, w))
end


"""
	sauterschwab_parameterized(integrand, method::SauterSchwabStrategy)

Compute interaction integrals using the quadrature introduced in [1].

Here, `integrand` is the pull-back of the integrand into the parametric domain
of the two triangles that define the integration domain.

The second argument 'strategy' is an object whose type is for triangles one of

	- `CommonFace`
	- `CommonEdge`
	- `CommonVertex`
	- `PositiveDistance`

and for quadrilaterals one of

	- `CommonFaceQuad`
	- `CommonEdgeQuad`
	- `CommonVertexQuad`

according to the configuration of the two patches defining the domain of integration.
The constructors of these classes take a single argument `acc` that defines
the number of quadrature points along each of the four axes of the final
rectangular (ξ,η) integration domain (see [1], Ch 5).

Note that here we use for a planar triangle the representation:

    x = x[3] + u*(x[1]-x[3]) + v*(x[2]-x[3])

with `u` ranging from 0 to 1 and `v` ranging from 0 to 1-`u`. This parameter
domain and representation is different from the one used in [1].

[1] Sauter. Schwwab, 'Boundary Element Methods', Springer Berlin Heidelberg, 2011
"""
function sauterschwab_parameterized(integrand, strategy::S) where {S <: SauterSchwabStrategy}

    qps = strategy.qps
    η1, w1 = qps[1]
    η2, w2 = qps[1]
    η3, w3 = qps[1]
    ξ, w4 = qps[1]
    acc = zero(w1 * w2 * w3 * w4 * strategy(integrand, η1, η2, η3, ξ))

    for (η1, w1) in qps
        for (η2, w2) in qps
            w12 = w1 * w2
            for (η3, w3) in qps
                w123 = w12 * w3
                for (ξ, w4) in qps
                    acc += w123 * w4 * strategy(integrand, η1, η2, η3, ξ)
                end
            end
        end
    end

    return acc
end
