""""
`Cochain(simplices::AbstractVector, values::AbstractVector)`

A type representing a cochain, i.e. a way to associate a value to each simplex. Simplices are expressed as an abstract vector `simplices` (generally represented as tuples). Values are stored in `values`. In practice, `values` will often be a sparse vector, in which case `findnz` returns the simplices where the cochian has non-zero value and the list of values on those simplices.

Use `keys` and `values` to access simplices and values respectively. The value at a specific simplex can be accessed with `getindex`.

## Examples

```julia
julia> c = Cochain([(1, 2), (2, 3), (1, 3)], [6, 9, -2])
Int64-valued 1-Cochain

julia> keys(c)
3-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (2, 3)
 (1, 3)

julia> values(c)
3-element Array{Int64,1}:
  6
  9
 -2

julia> c[(2, 3)]
9

julia> using SparseArrays

julia> c = Cochain([(1, 2), (2, 3), (1, 3), (1, 4)], sparse([6, 0, -2, 12]))
Sparse Int64-valued 1-Cochain. Non-zero values:
(1, 2) => 6
(1, 3) => -2
(1, 4) => 12
```
"""
struct Cochain{S<:AbstractVector, V<:AbstractVector}
    simplices::S
    values::V
end

Base.keys(c) = c.simplices
Base.values(c::Cochain) = c.values
Base.length(c::Cochain) = length(values(c))
Base.getindex(s::Cochain, t) = values(s)[searchsortedfirst(keys(s), t)]
Base.eltype(c::Type{Cochain{S, V}}) where {S, V} = eltype(V)
Base.eltype(c::Cochain) = eltype(typeof(c))

SparseArrays.issparse(c::Cochain) = issparse(values(c))
function SparseArrays.findnz(c::Cochain)
    issparse(c) || throw(ArgumentError("Argument of findnz must be sparse"))
    inds, vals = findnz(values(c))
    (keys(c)[inds], vals)
end

function Base.show(io::IO, s::Cochain)
    npts = nv(eltype(keys(s)))
    summary = "$(eltype(s))-valued $(npts-1)-Cochain"
    if !issparse(s)
        print(io, summary)
    else
        println(io, "Sparse $summary. Non-zero values:")
        limit = get(io, :limit, true)
        for (ind, (key, val)) in enumerate(zip(findnz(s)...))
            if ind >= 100 && limit
                println("...")
                break
            else
                show(io, key => val)
                println(io, "")
            end
        end
    end
end

function onecochain(::Type{T}, c::AbstractVector, i) where {T}
    v = sparsevec([i], one(T), length(c))
    Cochain(c, v)
end

@generated function boundaryvalue(cochain, simplex::T) where T
    N = nv(T)
    expr = :(cochain[face(simplex, 1)])
    for i in 2:N
        sign = isodd(i) ? (:+) : (:-) 
        expr = Expr(:call, sign, expr, :(cochain[face(simplex, $i)]))
    end
    expr
end
