
# Introduction

This package provides the Sauter-Schwab regularizing coordinate transformations [1] such that 4D integrals of the form
```math
\int_{\Gamma}\int_{\Gamma'}b_i(\bm{x})\,k(\bm{x},\bm{y})\, b_j(\bm{y})\,\mathrm{d}S(\bm{y})\,\mathrm{d}S(\bm{x})
```
with Cauchy-singular integral kernels ``k(\bm{x},\bm{y})`` can be integrated via numerical quadrature.
The integrals denote double surface integrals over 
- triangles (curved or flat) or 
- quadrilaterals (curved or flat) 
``\Gamma`` and ``\Gamma'`` in 3D Space. 
The functions ``b_i(\bm{x})`` and ``b_i(\bm{y})`` are assumed to be real valued and non-singular.

These kind of integrals occur in the area of boundary element methods (BEM) for solving elliptic partial differential equations. 
It can be interpreted as the interaction of the two basisfunctions ``b_i(\bm{x})`` and ``b_i(\bm{y})``, with respect to their domains ``\Gamma`` and ``\Gamma'``, which, for instance, correspond to the cells of a meshed surface.

!!! info
    The triangles or quadrilaterals must be either equal, have two vertices in common, have one vertex in common or do not touch at all. A partial overlap is forbidden.

    In the current implementation ``\Gamma`` and ``\Gamma'`` have to be both either triangles or quadrilatersls. However, mixed cases can be implemented, too.


## References

[1] Sauter S. Schwab C., "Boundary Element Methods (Springer Series in Computational Mathematics)", Chapter 5, Springer, 2010.
