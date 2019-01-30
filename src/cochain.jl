struct Cochain{S<:StructVector, V<:AbstractVector}
    simplices::S
    values::V
end

Base.keys(c) = c.simplices
Base.values(c::Cochain) = c.values
Base.length(c::Cochain) = length(values(c))
Base.getindex(s::Cochain, t::Tuple) = values(s)[searchsortedfirst(keys(s), t)]

Cochain(c::Cochain, v::AbstractVector) = Cochain(keys(c), v)

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
