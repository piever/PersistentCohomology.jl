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

persistentcocycles(::Type{T}, cplx, max_dim) where {T} =
    persistentcocycles(T, eltype(cplx[1]), cplx, max_dim)

persistentcocycles(::Type{T}, mat::AbstractSparseMatrix, max_dim) where {T} =
    persistentcocycles(T, vietorisrips(mat, max_dim+1), max_dim)
