# SauterSchwabQuadrature.jl  :-)


This package can be used to solve problems of following type:

```math
t(\vt{k},\vt{j}) = \frac{1}{ik} \int_{\Gamma} \int_{\Gamma'} \nabla \cdot \vt{k}(x) \nabla \cdot \vt{j}(y) \frac{e^{-ik|x-y|}}{4\pi|x-y|} dy dx - ik \int_{\Gamma} \int_{\Gamma'} \vt{k}(x) \cdot \vt{j}(y) \frac{e^{-ik|x-y|}}{4\pi|x-y|} dy dx
```

The above integral is a double area-integral over two triangles `Γ1` and `Γ2`. The integrand consists of two basisfunctions `b_i(x) ` and `b_j(y)` and the kernel `G(x,y)`.   

This kind of integral often appears in the area of Boundary Element Methods for solving elliptic partial differential equations and can be interpretated as the interaction of the two basisfunctions on their respective triangles.

As the solving algorithm works for a wide range of basisfunctions and kernels, all the conditions for the kernel, basisfunctions and the integration areas will be given:  

## Subtitle

Hello



## Index


```@index
```




```@docs
f(x)
```
