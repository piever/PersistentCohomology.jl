## Fallback

face(s::T, i) where {T} = error("face not implemented for type $T")

extend(s::T, i) where {T} = error("extend not implemented for type $T")

nv(::Type{T}) where {T} = error("nv not implemented for type $T")

## StaticArrays implementation

face(s::StaticVector, i) = deleteat(s, i)

extend(s::StaticVector, i) = push(s, i)

nv(::Type{<:StaticVector{N}}) where {N} = N

## Tuple implementation

function face(s::Tuple, i)
    Tuple(val for (ind, val) in enumerate(s) if ind != i)
end

extend(s::Tuple, i) = (s..., i)

nv(::Type{<:NTuple{N, Any}}) where {N} = N

