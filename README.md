# PersistentCohomology

[![Build Status](https://travis-ci.org/piever/PersistentCohomology.jl.svg?branch=master)](https://travis-ci.org/piever/PersistentCohomology.jl)
[![codecov.io](http://codecov.io/github/piever/PersistentCohomology.jl/coverage.svg?branch=master)](http://codecov.io/github/piever/PersistentCohomology.jl?branch=master)

Pure Julia package to compute persistent cohomology. In particular it implements the persistent cocycle algorithm from [Persistent Cohomology and Circular Coordinates - Vin de Silva, Dmitriy Morozov, Mikael Vejdemo-Johansson](https://link.springer.com/article/10.1007/s00454-011-9344-x). The algorithm can be used with coefficients on an arbitrary field: finite fields are available in Julia via the [GaloisFields](https://github.com/tkluck/GaloisFields.jl) package. 

This package provides a type to hold cochains (called `Cochain`) and two functions: `vietorisrips` to compute the filtered Vietoris Rips complex from a sparse distance matrix and `persistentcocycles` to compute persistent cohomolofy on such complex. See the docstrings of `Cochain`, `vietorisrips` and `persistentcocycles` for more details.

## Install

To install this package (and GaloisFields which you'll probably need as well):

```julia
pkg> add https://github.com/piever/PersistentCohomology.jl.git

pkg> add GaloisFields
```

## Examples

```julia
julia> using PersistentCohomology, SparseArrays, GaloisFields

julia> M = sparse(diagm(1 => fill(sqrt(2), 3), 2 => fill(2.0, 2), 3 => [sqrt(2)]))
4×4 SparseMatrixCSC{Float64,Int64} with 6 stored entries:
  [1, 2]  =  1.41421
  [1, 3]  =  2.0
  [2, 3]  =  1.41421
  [1, 4]  =  1.41421
  [2, 4]  =  2.0
  [3, 4]  =  1.41421

julia> cplx = vietorisrips(M, 2)
(Float64-valued 0-Cochain, Float64-valued 1-Cochain, Float64-valued 2-Cochain)

julia> const F3 = @GaloisFields 3

julia> cocycles, spans = persistentcocycles(F3, cplx, 1);

julia> cocycles[1]
4-element Array{Cochain{StructArrays.StructArray{Tuple{Int64},1,NamedTuple{(:x1,),Tuple{Array{Int64,1}}}},SparseVector{𝔽₃,Int64}},1}:
 Sparse 𝔽₃-valued 0-Cochain. Non-zero values:
(1,) => 1
(2,) => 1
(3,) => 1
(4,) => 1

 Sparse 𝔽₃-valued 0-Cochain. Non-zero values:
(2,) => 1

 Sparse 𝔽₃-valued 0-Cochain. Non-zero values:
(3,) => 1

 Sparse 𝔽₃-valued 0-Cochain. Non-zero values:
(4,) => 1


julia> spans[1]
4-element Array{IntervalSets.Interval{:closed,:open,Float64},1}:
 0.0..Inf (closed–open)
 0.0..1.4142135623730951 (closed–open)
 0.0..1.4142135623730951 (closed–open)
 0.0..1.4142135623730951 (closed–open)

julia> cocycles[2]
1-element Array{Cochain{StructArrays.StructArray{Tuple{Int64,Int64},1,NamedTuple{(:x1, :x2),Tuple{Array{Int64,1},Array{Int64,1}}}},SparseVector{𝔽₃,Int64}},1}:
 Sparse 𝔽₃-valued 1-Cochain. Non-zero values:
(4, 3) => 1


julia> spans[2]
1-element Array{IntervalSets.Interval{:closed,:open,Float64},1}:
 1.4142135623730951..2.0 (closed–open)
```
