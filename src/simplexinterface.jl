## Fallback

face(s::T, i) where {T} = error("face not implemented for type $T")

extend(s::T, i) where {T} = error("extend not implemented for type $T")

nv(::Type{T}) where {T} = error("nv not implemented for type $T")

## Tuple implementation

# From deleteat method in StaticArrays
@generated function face(vec::NTuple{N, Any}, index) where {N}
    newlen = N - 1
    exprs = [:(ifelse($i < index, vec[$i], vec[$i+1])) for i = 1:newlen]
    return quote
        Base.@_propagate_inbounds_meta
        @boundscheck if (index < 1 || index > $N)
            throw(BoundsError(vec, index))
        end
        @inbounds return tuple($(exprs...))
    end
end

extend(s::Tuple, i) = (s..., i)

nv(::Type{<:NTuple{N, Any}}) where {N} = N

