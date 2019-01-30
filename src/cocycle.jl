const Span{U} = Interval{:closed, :open, U}

function initialize_cocycles(::Type{T}, ::Type{U}, cplx) where {T, U}
    instances = map(t -> Cochain(keys(t), spzeros(T, length(t))), cplx)
    map(instances) do inst
        typ = NamedTuple{(:cocycle, :span), Tuple{typeof(inst), Span{U}}}
        StructVector{typ}(undef, 0)
    end
end

function update_cocycles!(::Type{T}, ::Type{U}, simplices,
                          c_small, c_big, indices, weight) where {T, U}
    for index in indices
        simplex = simplices[index]
        m = map(c_small) do (coc, interval)
            weight in interval ? boundaryvalue(coc, simplex) : zero(U)
        end
        s = findall(!iszero, m)
        if !isempty(s)
            li = last(s)
            lc, ic = c_small[li]
            for ii in s[1:end-1]
                ratio = m[ii] / m[li]
                vals = values(c_small[ii].cocycle)
                vals .-= ratio .* values(lc)
            end
            t0, _ = endpoints(ic)
            c_small.span[li] = Span{U}(t0, weight)
        elseif c_big !== nothing
            cocycle = onecochain(T, simplices, index)
            span = Span{U}(weight, U(Inf))
            push!(c_big, (cocycle = cocycle, span = span))
        end
    end
end

function persistentcocycles(::Type{T}, ::Type{U}, cplx, max_dim) where {T, U}
    cocycles_extended = initialize_cocycles(T, U, cplx)
    cocycles = cocycles_extended[1:min(length(cocycles_extended), max_dim+1)] 
    sv = simplexiterator(cplx)
    for (weight, (npts, indices)) in sv
        npts > length(cocycles)+1 && continue
        c_small = npts == 1 ? () : cocycles[npts-1]
        c_big = npts <= length(cocycles) ? cocycles[npts] : nothing
        update_cocycles!(T, U, keys(cplx[npts]), c_small, c_big, indices, weight)
    end
    fc = map(cocycles) do c
        filter(t -> !isempty(t.span), c)
    end
    (map(t -> t.cocycle, fc), map(t -> t.span, fc))
end

"""
`persistentcocycles(T, cplx, max_dim = 1)`

Implementation of the persistent cocycle algorithm from Persistent Cohomology and Circular Coordinates - Vin de Silva, Dmitriy Morozov, Mikael Vejdemo-Johansson (https://link.springer.com/article/10.1007/s00454-011-9344-x). 

Type `T` must represent a field. Finite fields are available in Julia via the GaloisFields package. `cplx` is a filtered complex (where simplices are assumed to be sorted): see `vietorisrips` for a description of such complex and a way to build it from a sparse distance matrix.

Return two tuples: cocycles (in turn represented as a tuple from 0-cocycles to `max_dim`-cocycles) and spans (again as a tuple) representing the interval during which each cocycle lives.

## Examples

```julia
julia> using SparseArrays

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

julia> using GaloisFields

julia> cocycles, spans = persistentcocycles(@GaloisField(3), cplx, 1);

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
"""
persistentcocycles(::Type{T}, cplx, max_dim = 1) where {T} =
    persistentcocycles(T, eltype(cplx[1]), cplx, max_dim)

persistentcocycles(::Type{T}, mat::AbstractSparseMatrix, max_dim = 1) where {T} =
    persistentcocycles(T, vietorisrips(mat, max_dim+1), max_dim)
