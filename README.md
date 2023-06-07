# LinearAlgebra.jl

Implements basic functionality for linear algebra: Matrix operations, Gaussian Elimination, LU Decomposition

## Usage
```julia
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.7.2 (2022-02-06)
 _/ |\__'_|_|_|\__'_|  |  HEAD/bf53498635 (fork: 461 commits, 536 days)
|__/                   |

pkg> add https://github.com/jin-park/linearalgebra.jl.git
pkg> using LinearAlgebra
```

Press "]" to go into pkg mode in the Julia REPL

## Relevant functions

```julia
gaussian_elimination(A, b)
LU(A::AbstractMatrix)
sys_of_linear_eqs(A, b)
```

## Note

This package is unoptimized for performance and probably too slow for any real world usage. Publishing on github just for fun ;) 
