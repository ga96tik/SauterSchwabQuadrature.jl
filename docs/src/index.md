# SauterSchwabQuadrature.jl  THIS IS JUST A TEST


This package can be used to solve problems of following type:

```math
```
The above integral is a double area-integral over two triangles `Γ1` and `Γ2`. The integrand consists of two basisfunctions `b_i(x) ` and `b_j(y)` and the kernel `G(x,y)`.   

This kind of integral often appears in the area of Boundary Element Methods for solving elliptic partial differential equations and can be interpretated as the interaction of the two basisfunctions on their respective triangles.

As the solving algorithm works for a wide range of basisfunctions and kernels, all the conditions for the kernel, basisfunctions and the integration areas will be given:  

## Subtitle

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```



## Index


```@index
```




```@docs
f(x)
```
