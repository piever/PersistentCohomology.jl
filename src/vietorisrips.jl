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
    push!(cocycle.simplices, sim)
    push!(cocycle.values, sim_weight)

    (length(sim) > dim_max) && return

    for (n, w) in weighted_neighbors
        ptr = g.colptr[n]:g.colptr[n+1]-1
        new_weighted_neighbors = neighbors_intersection(g.rowval[ptr], g.nzval[ptr], weighted_neighbors)
        add_cofaces!(v, g, (n, sim...), max(w, sim_weight), new_weighted_neighbors, dim_max)
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
