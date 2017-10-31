# SauterSchwabQuadrature.jl

This package can be used to solve problems of following type:

```math
```
The above integral is a double area-integral over two triangles `Γ1` and `Γ2`. The integrand consists of two basisfunctions `b_i(x) ` and `b_j(y)` and the kernel `G(x,y)`.   

This kind of integral often appears in the area of Boundary Element Methods for solving elliptic partial differential equations and can be interpretated as the interaction of the two basisfunctions on their respective triangles.

As the solving algorithm works for a wide range of basisfunctions and kernels, all the conditions for the kernel, basisfunctions and the integration areas will be given:   


## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs help` - Print this help message.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
