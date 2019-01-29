const Span{U} = Interval{:closed, :open, U}

function initialize_cocycles(::Type{T}, ::Type{U}, cplx) where {T, U}
    instances = map(t -> Cochain(t, spzeros(T, length(t))), cplx)
    map(instances) do inst
        typ = NamedTuple{(:cocycle, :span), Tuple{typeof(inst), Span{U}}}
        StructVector{typ}(undef, 0)
    end
end

function generate_remove_cocycles!(::Type{T}, ::Type{U}, simplices, c_small, c_big, row) where {T, U}
    for index in row.indices
        simplex = simplices[index]
        m = map(c_small) do (coc, interval)
            row.weight in interval ? boundaryvalue(coc, simplex) : zero(U)
        end
        s = findall(!iszero, m)
        if !isempty(s)
            li = last(s)
            lc, ic = c_small[li]
            for ii in s[1:end-1]
                ratio = m[ii] / m[li]
                c_small[ii].cocycle.values .-= ratio .* lc.values
            end
            t0, _ = endpoints(ic)
            c_small.span[li] = Span{U}(t0, row.weight)
        elseif c_big !== nothing
            cocycle = onecochain(T, simplices, index)
            span = Span{U}(row.weight, U(Inf))
            push!(c_big, (cocycle = cocycle, span = span))
        end
    end
end

function persistent_cocycles(::Type{T}, ::Type{U}, cplx, max_dim) where {T, U}
    cocycles_extended = initialize_cocycles(T, U, cplx)
    cocycles = cocycles_extended[1:min(length(cocycles_extended), max_dim+1)] 
    sv = simplexiterator(cplx)
    for (weight, (npts, indices)) in sv
        npts > length(cocycles)+1 && continue
        row = (indices = indices, npts = npts, weight = weight)
        c_small = npts == 1 ? () : cocycles[npts-1]
        c_big = npts <= length(cocycles) ? cocycles[npts] : nothing
        generate_remove_cocycles!(T, U, cplx[npts].simplices, c_small, c_big, row)
    end
    map(cocycles) do c
        filter(t -> !isempty(t.span), c)
    end
end

persistent_cocycles(::Type{T}, cplx, max_dim) where {T} =
    persistent_cocycles(T, eltype(cplx[1].values), cplx, max_dim)

