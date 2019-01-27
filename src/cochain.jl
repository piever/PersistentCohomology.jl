struct Cochain{S<:StructVector, V<:AbstractVector}
    simplices::S
    values::V
end

Base.length(c::Cochain) = length(c.values)
Base.getindex(s::Cochain, t::Tuple) = s.values[searchsortedfirst(s.simplices, t)]

# From IndexedTables to get similar functionality as NDSparse
Base.broadcast(f::Function, s::Cochain, t::Cochain) = Cochain(s.simplices, broadcast(f, s.values, t.values))
Base.broadcast(f::Function, x::Cochain) = Cochain(x.simplices, broadcast(f, x.values))
Base.broadcast(f::Function, x::Cochain, y) = Cochain(x.simplices, broadcast(f, x.values, y))
Base.broadcast(f::Function, y, x::Cochain) = Cochain(x.simplices, broadcast(f, y, x.values))

Broadcast.broadcasted(f::Function, A::Cochain) = broadcast(f, A)
Broadcast.broadcasted(f::Function, A::Cochain, B::Cochain) = broadcast(f, A, B)
Broadcast.broadcasted(f::Function, A, B::Cochain) = broadcast(f, A, B)
Broadcast.broadcasted(f::Function, A::Cochain, B) = broadcast(f, A, B)

function weightcochain(cells, weights)
    s = StructArray(Tuple(a[:values]) for a in cells)
    p = sortperm(s)
    sorted = s[p]
    Cochain(sorted, weights[p])
end

function weightcochains(cplx, w)
    Dict(key => weightcochain(cells(cplx, key), val) for (key, val) in w)
end

function face(s, i)
    Tuple(val for (ind, val) in enumerate(s) if ind != i)
end

boundaryvalue(t::NTuple{N, Any}, ::Missing) where {N} = missing

@generated function boundaryvalue(t::NTuple{N, Any}, cells) where N
    expr = :(cells[face(t, 1)])
    for i in 2:N
        sign = isodd(i) ? (:+) : (:-) 
        expr = Expr(:call, sign, expr, :(cells[face(t, $i)]))
    end
    expr
end

function simplexlist(cochains)
    sv = vcat((StructArray(weight = val.values, dim = fill(key, length(val)), index = 1:length(val)) for (key, val) in cochains)...)
    sort!(sv)
    sv
end
