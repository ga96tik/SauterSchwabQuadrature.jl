# SauterSchwabQuadrature.jl


This package can be used to solve problems of following type:

```math
\int_{\Gamma}\int_{\Gamma'}b_i(\textbf{x})\,k(\textbf{x},\textbf{y})\, b_j(\textbf{y})\;da_\textbf{y}\,da_\textbf{x}
```

The above expression is a double area-integral over two flat triangles $\Gamma$ and $\Gamma'$ in 3D Space. The integrand consists of two basisfunctions $b_i(\textbf{x})$ and $b_i(\textbf{y})$ and the kernel $k(\textbf{x},\textbf{y})$.   

This kind of integral often appears in the area of Boundary Element Method for solving elliptic partial differential equations, and can be interpretated as the interaction of the two basisfunctions with respect to their triangles. For this reason in this package the two triangles are called test- and sourcecell, and the same goes for the two basefunctions, they are called test- and sourcefunction. The triangles correspond to the cells of a meshed surface.

As the solving algorithm works for a wide range of basisfunctions and kernels, all the requirements for the kernel, basisfunctions and the integration areas will be given:

1.Requirements for the triangles:
* The triangles must be flat and are created by three vertices
* The triangles must be either the same, have two vertices in common, have one vertex in common or do not touch at all; A partial overlap is forbidden

2.Requirements for the basisfunctions:
* The basisfunctions must be real and non singular on their respective triangles
* The basisfunctions map vectors on scalars

3.The kernel must be Cauchy-Singular

According to item 1, four different constellations of the two triangles are possible:
* Equal triangles $\to$ Common Face
* Two vertices in common $\to$ Common Edge
* One vertex in commmon $\to$ Common Vertex
* Both triangles do not touch at all $\to$ Positive Distance

As each of those four constellations has its own integration method (because of a possible singularity in the kernel), every single case will be presented. However, in all four cases the called function will look like:

`sauterschwabintegral(sourcechart, testchart, integrand, constellation)`.

`sourcechart` and `testchart` are the mappings from a reference triangle (parametrisation) to the real triangles in space, they contain the information of the positions of both triangles. `integrand` is the integrand and `constellation` is a object of a particular type which conveys the function the method of integration and additional information of how accurate the integration shall be done.

As this package depends on the package 'CompScienceMeshes', it is useful to be familiar with that package, at least with the functions `simplex()` and `point()`.

This documentation does not derive the integration rules, and how the integration is done, it only shows how to handle this package. If the user wants to know more about how this package operates, he has to go inside the src folder and look up for the book in the README file.






## Common Face

 $\Gamma$ and $\Gamma'$ are equal, hence `sourcechart` and `testchart` are equal as well. The two charts can be created by

 `testchart = sourcechart = simplex(P1,P2,P3)
 `

 , where `P1`, `P2` and `P3` are the vertices of that particular triangle. Note that both charts must be equal, that means that the first argument of both charts must be euqal, the second argument of both charts must be equal, and the last argument of both charts must be equal. For instance `sourcechart = simplex(P1,P2,P3)` and  `testchart = simplex(P1,P3,P2)` will not work because of the different orders of the input arguments.

 The `integrand` must be defined as a function and its function name is the input argument.

 The last argument can be created by

`cf = CommonFace(x)`.

`cf` is an object of type `CommonFace()`, x is an integer which stands for the degree of accuracy. The larger x is, the more accurate the evaluation will be.

An example of this case can be found and run in the example folder.






## Common Edge

$\Gamma$ and $\Gamma'$ are now different, hence `sourcechart` and `testchart` are different as well. The two charts have to be created in the following manner:

`testchart = simplex(P1,P2,P3); sourcechart = simplex(P1,P4,P3)`.

Again the order of the input arguments must be taken into account: The first argument of both charts must be equal, and the last argument of both charts must be equal. Consequently the first and the last argument are the vertices which both triangles have in common.

The `integrand` must be defined as a function and its function name is the input argument.

The last argument can be created by

`ce = CommonEdge(x)`.

`ce` is an object of type `CommonEdge()`, x is an integer which stands for the degree of accuracy. The larger x is, the more accurate the evaluation will be.

An example of this case can be found and run in the example folder.






## Common Vertex

The two triangles and charts are again different. The two charts have to be created in the following manner:

`sourcechart = simplex(P1,P2,P3); testchart = simplex(P1,P4,P5)`.

Again the order of the input arguments must be taken into account: The first argument of both charts must be equal, the order of `P2` and `P3` with respect to `sourcechart` and the order of `P4` and `P5` with respect to `testchart` does not matter.  Consequently, the first argument is the vertex both triangles have in common.

The `integrand` must be defined as a function and its function name is the input argument.

The last argument can be created by

`cv = CommonVertex(x)`.

`cv` is an object of type `CommonVertex()`, x is an integer which stands for the degree of accuracy. The larger x is, the more accurate the evaluation will be.

An example of this case can be found and run in the example folder.






## Positive Distance

As the triangles do not touch at all, the integration becomes a simple quadrature. Therefore the order of arguments for the two `simplexfunction()`'s no longer matter.

The `integrand` must be defined as a function and its function name is the input argument.

The last argument can be created by

`pd = PositiveDistance(x)`.

`pd` is an object of type `PositiveDistance()`, x is an integer which stands for the degree of accuracy. The larger x is, the more accurate the evaluation will be.

An example of this case can be found and run in the example folder.
