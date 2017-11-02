
<a id='SauterSchwabQuadrature.jl-THIS-IS-JUST-A-TEST-1'></a>

# SauterSchwabQuadrature.jl  THIS IS JUST A TEST


This package can be used to solve problems of following type:


$$

$$


The above integral is a double area-integral over two triangles `Γ1` and `Γ2`. The integrand consists of two basisfunctions `b_i(x)` and `b_j(y)` and the kernel `G(x,y)`.   


This kind of integral often appears in the area of Boundary Element Methods for solving elliptic partial differential equations and can be interpretated as the interaction of the two basisfunctions on their respective triangles.


As the solving algorithm works for a wide range of basisfunctions and kernels, all the conditions for the kernel, basisfunctions and the integration areas will be given:  


<a id='Subtitle-1'></a>

## Subtitle


$$
\frac{n!}{k!(n - k)!} = \binom{n}{k}
$$


<a id='Index-1'></a>

## Index

- [`SauterSchwabQuadrature.f`](index.md#SauterSchwabQuadrature.f-Tuple{Any})

<a id='SauterSchwabQuadrature.f-Tuple{Any}' href='#SauterSchwabQuadrature.f-Tuple{Any}'>#</a>
**`SauterSchwabQuadrature.f`** &mdash; *Method*.



f(x) is a super duper cool function


<a target='_blank' href='https://github.com/ga96tik/SauterSchwabQuadrature.jl/tree/d1a0975fe1f96e96777df19c4d61f7fe4d849c6e/src/f.jl#L4-L6' class='documenter-source'>source</a><br>

