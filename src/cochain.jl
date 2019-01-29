struct Cochain{S<:StructVector, V<:AbstractVector}
    simplices::S
    values::V
end

Base.length(c::Cochain) = length(c.values)
Base.getindex(s::Cochain, t::Tuple) = s.values[searchsortedfirst(s.simplices, t)]

Cochain(c::Cochain, v::AbstractVector) = Cochain(c.simplices, v)

function onecochain(::Type{T}, c::StructVector, i) where {T}
    null = spzeros(T, length(c))
    null[i] = one(T)
    Cochain(c, null)
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
