# SauterSchwabQuadrature.jl


This package can be used to solve problems of following type:

```math
\int_{\Gamma}\int_{\Gamma'}b_i(\textbf{x})\,k(\textbf{x},\textbf{y})\, b_j(\textbf{y})\;da_\textbf{y}\,da_\textbf{x}
```

The above expression is a double area-integral over two triangles (curved or flat) $\Gamma$ and $\Gamma'$ in 3D Space. The integrand consists of two basisfunctions $b_i(\textbf{x})$ and $b_i(\textbf{y})$ and the kernel $k(\textbf{x},\textbf{y})$.   

This kind of integral appears in the area of Boundary Element Methods for solving elliptic partial differential equations, and can be interpretated as the interaction of the two basisfunctions with respect to their triangles. For this reason in this package the two triangles are called test- and sourcecell, and the same goes for the two basefunctions; they are called test- and sourcefunction. The triangles correspond to the cells of a meshed surface.

As the solving algorithm works for a wide range of basisfunctions and kernels, all the requirements for the kernel, basisfunctions and the integration areas will be given:

1.Requirements for the triangles:
* The triangles are created by three vertices
* The triangles must be either the same, have two vertices in common, have one vertex in common or do not touch at all; A partial overlap is forbidden

2.Requirements for the basisfunctions:
* The basisfunctions must be real and non singular on their respective triangles
* The basisfunctions map vectors on scalars

3.The kernel must be Cauchy-Singular

Depending on the input data, two different implementations of the integral are contained in this package. The first one is very convenient to handle, and does not need a parametrisation; but it works only for flat triangles, and the user has to be familiar with the functions `simplex()` and `point()` of the package CompScienceMeshes. The second implementation contains only the integration rules; so the user has to build the parameterization himself, but therefore it also works for curved triangles.

Both implementations are called by a function, which looks like:  

`function(sourcechart, testchart, integrand, information)`

`sourcechart` and `testchart` are the mappings from a reference triangle (parameterization) to the real triangles in space; they contain the information of the positions of both triangles. `integrand` is the integrand, and the last argument contains information about how accurate the integration shall be done, and the type of integration.

On the pages 'Non-Parameterized' and 'Parameterized', the user will find more information about the two implementations, and how to operate those.

Before using this package, CompScienceMeshes has to be installed; even if the user does not need it (second implementation), this package depends internally on CompScienceMeshes.  

This documentation does not derive the integration rules, and how the integration is done, it only shows how to handle this package. If the user wants to know more about how this package operates, he has to go inside the src folder and look up for the book quoted in the README file.
