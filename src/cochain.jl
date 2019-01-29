struct Cochain{S<:StructVector, V<:AbstractVector}
    simplices::S
    values::V
end

Base.length(c::Cochain) = length(c.values)
Base.getindex(s::Cochain, t::Tuple) = s.values[searchsortedfirst(s.simplices, t)]

Cochain(c::Cochain, v::AbstractVector) = Cochain(c.simplices, v)

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

function simplexlist(cochains)
    sv = vcat((StructArray(weight = val.values, npts = fill(ind, length(val)), index = 1:length(val)) for (ind, val) in enumerate(cochains))...)
    sort!(sv)
    sv
end
