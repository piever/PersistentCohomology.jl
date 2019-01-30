struct Cochain{S<:StructVector, V<:AbstractVector}
    simplices::S
    values::V
end

Base.keys(c) = c.simplices
Base.values(c::Cochain) = c.values
Base.length(c::Cochain) = length(values(c))
Base.getindex(s::Cochain, t::Tuple) = values(s)[searchsortedfirst(keys(s), t)]
Base.eltype(c::Type{Cochain{S, V}}) where {S, V} = eltype(V)
Base.eltype(c::Cochain) = eltype(typeof(c))

SparseArrays.issparse(c::Cochain) = issparse(values(c))
function SparseArrays.findnz(c::Cochain)
    issparse(c) || throw(ArgumentError("Argument of findnz must be sparse"))
    inds, vals = findnz(values(c))
    (keys(c)[inds], vals)
end

function Base.show(io::IO, s::Cochain)
    npts = length(propertynames(keys(s)))
    summary = "$(eltype(s))-valued Cochain of dimension $(npts-1)"
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

function onecochain(::Type{T}, c::StructVector, i) where {T}
    v = sparsevec([i], one(T), length(c))
    Cochain(c, v)
end

function face(s, i)
    Tuple(val for (ind, val) in enumerate(s) if ind != i)
end

@generated function boundaryvalue(cochain, simplex::NTuple{N, Any}) where N
    expr = :(cochain[face(simplex, 1)])
    for i in 2:N
        sign = isodd(i) ? (:+) : (:-) 
        expr = Expr(:call, sign, expr, :(cochain[face(simplex, $i)]))
    end
    expr
end
