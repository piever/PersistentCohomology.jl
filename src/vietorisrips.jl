"""
`vietorisrips(g::AbstractSparseMatrix{U, T}, dim_max = 2) where {U, T}`

Return a filtered Vietoris Rips complex (up to dimension `dim_max`) from a sparse distance matrix `g`. Only the upper triangular part of `g` is used. Zero values of the matrix are ignored: to obtain the rips filtration up to a threshold `t`, simply remove all entries of `g` greater than `t`. Discretizing values of entries in `g` may also improve performance.

The filtered complex is a tuple (where the first element corresponds to 0-simplices, the second element to 1-simplices and so on). The filtration on n-simplices is represented (slightly improperly) as a `U`-valued n-`Cochain`.

## Examples

```julia
julia> using SparseArrays

julia> M = sparse([1, 2, 1], [2, 3, 3], [0.3, 0.2, 0.4], 3, 3)
3×3 SparseMatrixCSC{Float64,Int64} with 3 stored entries:
  [1, 2]  =  0.3
  [1, 3]  =  0.4
  [2, 3]  =  0.2

julia> cplx = vietorisrips(M, 2)
(Float64-valued 0-Cochain, Float64-valued 1-Cochain, Float64-valued 2-Cochain)

julia> keys(cplx[2])
3-element StructArrays.StructArray{Tuple{Int64,Int64},1,NamedTuple{(:x1, :x2),Tuple{Array{Int64,1},Array{Int64,1}}}}:
 (2, 1)
 (3, 1)
 (3, 2)

julia> values(cplx[2])
3-element Array{Float64,1}:
 0.3
 0.4
 0.2
```
"""
function vietorisrips(g::AbstractSparseMatrix, dim_max = 2)
    lt = sparse_uppertriangular(g)
    vietorisrips_uppertriangular(lt, dim_max)
end

function sparse_uppertriangular(g::AbstractSparseMatrix)
    size1, size2 = size(g)
    @assert size1 == size2
    row, col, vals = findnz(g)
    mask = row .< col
    sparse(row[mask], col[mask], vals[mask], size1, size2)
end 

function vietorisrips_uppertriangular(g::SparseMatrixCSC{U, T}, dim_max = 2) where {U, T}
    v = ntuple(dim_max+1) do npts
        Cochain(StructVector{NTuple{npts, T}}(undef, 0), U[])
    end
    compute_simplices!(v, g, dim_max)
    return v
end

function compute_simplices!(v, g::SparseMatrixCSC{U, T}, dim_max) where {U, T}
    us = one(T):T(size(g, 1))
    for u in us
        ptr = g.colptr[u]:g.colptr[u+1]-1
        nu = view(g.rowval, ptr)
        ws = view(g.nzval, ptr)
        add_cofaces!(v, g, (u,), zero(U), StructVector((nu, ws)), dim_max)
    end
end

function add_cofaces!(v, g, sim, sim_weight, weighted_neighbors, dim_max)
    cocycle = v[length(sim)]
    push!(keys(cocycle), sim)
    push!(values(cocycle), sim_weight)

    (length(sim) > dim_max) && return

    for (n, w) in weighted_neighbors
        ptr = g.colptr[n]:g.colptr[n+1]-1
        new_weighted_neighbors = neighbors_intersection(g.rowval[ptr], g.nzval[ptr], weighted_neighbors)
        add_cofaces!(v, g, (sim..., n), max(w, sim_weight), new_weighted_neighbors, dim_max)
    end
end

function neighbors_intersection(neighbors, weights, weighted_neighbors)
    new_weighted_neighbors = similar(weighted_neighbors, 0)
    for (nb, nb_weight) in weighted_neighbors
        idx = searchsortedfirst(neighbors, nb)
        idx in axes(neighbors, 1) || continue
        neighbors[idx] == nb || continue
        push!(new_weighted_neighbors, (nb, max(nb_weight, weights[idx])))
    end
    new_weighted_neighbors
end
