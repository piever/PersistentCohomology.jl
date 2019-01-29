###
# Courtesy of @shashi

struct SortMergeItr{S, T, F<:Function}
    a::S
    b::T
    lt::F
end

SortMergeItr(a, b) = SortMergeItr(a, b, isless)

import Base: IsInfinite, SizeUnknown, HasLength, HasShape,
             IteratorSize, IteratorEltype, HasEltype, EltypeUnknown

function IteratorSize(::Type{SortMergeItr{S,T, F}}) where {S,T, F}
    iS = IteratorSize(S)
    iT = IteratorSize(T)
    if all(isa.((iS, iT), Union{HasLength, HasShape}))
        return HasLength()
    elseif IsInfinite() in (iS, iT)
        return IsInfinite()
    elseif SizeUnknown() in (iS, iT)
        return SizeUnknown()
    end
end
Base.length(itr::SortMergeItr) = length(itr.a) + length(itr.b)
Base.eltype(::Type{SortMergeItr{S,T, F}}) where {S,T, F} = Union{eltype(S), eltype(T)}
function IteratorEltype(::Type{SortMergeItr{S,T, F}}) where {S,T,F}
    if IteratorEltype(S) == HasEltype() && IteratorEltype(T) == HasEltype()
        HasEltype()
    else
        EltypeUnknown()
    end
end

function _pick(a, b, lt)
    if a === nothing && b === nothing
        return nothing
    elseif a === nothing
        return b[1], (2, a, b)
    elseif b === nothing
        return a[1], (1, a, b)
    else
        lt(a[1], b[1]) ? (a[1], (1, a, b)) : (b[1], (2, a, b))
    end
end

function Base.iterate(s::SortMergeItr)
    a = iterate(s.a)
    b = iterate(s.b)
    _pick(a, b, s.lt)
end

function Base.iterate(s::SortMergeItr, (select, a, b,))
    if select == 1
        a = iterate(s.a, a[2])
    else
        b = iterate(s.b, b[2])
    end

    _pick(a, b, s.lt)
end

###

function hassmallerkey(a, b)
    a_key, _ = a
    b_key, _ = b
    isless(a_key, b_key)
end

function simplexiterator(cplx)
    vals = map(t -> t.values, cplx)
    itrs = map(finduniquesorted, vals)
    aug_itrs = map(enumerate(itrs)) do (n, itr)
        (i[1] => (n, i[2]) for i in itr)
    end

    foldl((a, b) -> SortMergeItr(b, a, hassmallerkey), aug_itrs)
end

