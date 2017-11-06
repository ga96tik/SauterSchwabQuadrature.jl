# SauterSchwabQuadrature.jl


This package can be used to solve problems of following type:

```math
\int_{\Gamma}\int_{\Gamma'}b_i(\textbf{x})\,k(\textbf{x},\textbf{y})\, b_j(\textbf{y})\;da_\textbf{y}\,da_\textbf{x}
```

The above expression is a double area-integral over two flat triangles $\Gamma$ and $\Gamma'$ in 3D Space. The integrand consists of two basisfunctions $b_i(\textbf{x})$ and $b_i(\textbf{y})$ and the kernel $k(\textbf{x},\textbf{y})$.   

This kind of integral often appears in the area of Boundary Element Methods for solving elliptic partial differential equations, and can be interpretated as the interaction of the two basisfunctions on their respective triangles. Therefore in this package the two triangles are called test- and sourcecell, and the same goes for the two basefunctions, they are called test- and sourcefunction. The triangles correspond to the cells of a meshed surface.

As the solving algorithm works for a wide range of basisfunctions and kernels, all the requirements for the kernel, basisfunctions and the integration areas will be given:

1. Requiremts for the triangles:
* The triangles must be flat and are created by three vertices
* The triangles must be either the same, or have two vertices in common, or have one vertex in common or they do not touch at all; A partial overlap is forbidden

2. Requirements for the basisfunctions:
* The basisfunctions must be real and non singular on their respective triangles
* The basisfunctions map vectors on scalars

3. The kernel must be Cauchy-Singular

According to item 1, four different constellations of the two triangles are possible:
*  Equal triangles $\to$ Common Face
* Two vertices in common $\to$ Common Edge
* One vertex in commmon $\to$ Common Vertex
* Both triangles do not touch at all $\to$ Positive Distance

As each of those four constellations has its own integration method (because of a possible singularity in the kernel), every single method will be investigated in this documentation. However, in all four cases the called function will look like:
```
sauterschwabintegral(sourcechart, testchart, kernel, method)
```
`Sourcechart` and `Testchart` are the mappings from a reference triangle to the real triangles in space, they contain the information of the positions of both triangles. `kernel` is the integrand and `method` is a object of a particular type which conveys the function the method of integration and additionally it contains information of how accurate the integration shall be done.

As this package depends on CompScienceMeshes, it is useful to know about that package, at least about the function `simplex`.




## Subtitle

Hello



## Index


```@index
```




```@docs
f(x)
```
