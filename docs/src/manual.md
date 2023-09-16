
# Manual

Fundamentally, one function is provided:
```@julia
sauterschwab_parameterized(integrand, strategy)
```

The first argument `integrand` is the parameterized integrand ``k'(\chi_\tau(u,v), \chi_t(u',v'))``. 
That is, it takes as argument two tuples: 
```@julia
integrand((u,v), (u',v'))
```

The second argument `strategy` specifies the reparametrization and is one of (the by this package provided) structs:

For triangles
- `CommonFace`
- `CommonEdge`
- `CommonVertex`
- `PositiveDistance`


For quadrilaterals
- `CommonFaceQuad`
- `CommonEdgeQuad`
- `CommonVertexQuad`

Each such struct takes one argument specifying the quadrature rule, e.g.,
```@julia
strategy = CommonEdge(qrule)
```
where `qrule` is a vector of `(point, weight)` tuples for a quadrature on the domain ``[0,1]``.

!!! tip
    We recommend the [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl/tree/master) package.
    For a Gauss-Legendre quadrature a method is provided that maps to the ``[0,1]`` domain:
    ```@julia
    order = 10
    qrule = SauterSchwabQuadrature._legendre(order, 0, 1)
    ```