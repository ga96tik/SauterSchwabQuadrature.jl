
# Details

The integrals of the form
```math
\int_{\Gamma}\int_{\Gamma'}b_i(\bm{x})\,k(\bm{x},\bm{y})\, b_j(\bm{y})\,\mathrm{d}S(\bm{y})\,\mathrm{d}S(\bm{x})
```
are solved based on relative coordinates.

## Relative Coordinates

It is assumed that by suitable parameter transforms
```math
\chi_\tau: \hat\tau \mapsto \tau
```
with ``\hat\tau = (u,v) \in [0,1]^2`` and
```math
\chi_t: \hat t \mapsto t
```
with ``\hat t = (u',v') \in [0,1]^2`` from the reference element (triangle or square) to the actual ones, the integral is brought into the form
```math
\iint_{\hat \tau}\iint_{\hat t} k'(\chi_\tau(u,v), \chi_t(u',v'))  \, \mathrm{d}u \mathrm{d}v \mathrm{d}u' \mathrm{d}u' 
```
where ``k'(\bm{x},\bm{y})`` contains  ``k``, ``b_i``, ``b_j``, and the Jacobi determinant of the parametrization.

!!! tip
    The parametrizations for triangles and quadrilaterals are commonly employed for the integration without regularizations, as well.
    Hence, they are often already available.

For the regularizing parametertransform according to [1] four cases are distinguished.



### Common Face Case

``\Gamma`` and ``\Gamma'`` are equal, and both parameterizations must be equal, that is, ``\chi_t = \chi_\tau``.

```@raw html
<div align="center">
<img src="../assets/CommonFace.jpg" width="600"/>
</div>
<br/>
```


### Common Edge Case

``\Gamma`` and ``\Gamma'`` have an edge in common, and both parameterizations must fulfill the condition ``\chi_t(s,0) = \chi_\tau(s,0)``. For example, this condition could be met if the points ``(u\in[0,1];0)`` and ``(u'\in[0,1];0)`` are mapped on the same point on the common edge.

```@raw html
<div align="center">
<img src="../assets/CommonEdge.jpg" width="600"/>
</div>
<br/>
```


### Common Vertex Case

``\Gamma`` and ``\Gamma'`` have one vertex in common, and both parameterizations must fulfill the condition ``\chi_t(0,0) = \chi_\tau(0,0)``. 
This condition means, that the origin of both reference triangles is mapped on the common vertex.

```@raw html
<div align="center">
<img src="../assets/CommonVertex.jpg" width="600"/>
</div>
<br/>
```


### Positive Distance Case

The two triangles do not touch at all, and both parameterizations only need to map from the reference triangle onto the real triangle.

```@raw html
<div align="center">
<img src="../assets/PositiveDistance.jpg" width="600"/>
</div>
<br/>
```

