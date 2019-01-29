const Span{U} = Interval{:closed, :open, U}

function initialize_cocycles(::Type{T}, ::Type{U}, cplx) where {T, U}
    instances = map(t -> Cochain(t, spzeros(T, length(t))), cplx)
    map(instances) do inst
        typ = NamedTuple{(:cocycle, :span), Tuple{typeof(inst), Span{U}}}
        StructVector{typ}(undef, 0)
    end
end

function onecochain(::Type{T}, c::Cochain, i) where {T}
    null = spzeros(T, length(c))
    null[i] = one(T)
    Cochain(c, null)
end

function persistent_cocycles(::Type{T}, ::Type{U}, cplx, sv, max_dim) where {T, U}
    cocycles_extended = initialize_cocycles(T, U, cplx)
    cocycles = cocycles_extended[1:min(length(cocycles_extended), max_dim+1)] 
    for row in sv
        row.npts > length(cocycles)+1 && continue
        if row.npts == 1
            push!(cocycles[1], (cocycle = onecochain(T, cplx[1], row.index), span = Span{U}(zero(U), U(Inf))))
        else
            c_small = cocycles[row.npts-1]
            simplex = cplx[row.npts].simplices[row.index]
            m = map(c_small) do (coc, interval)
                row.weight in interval ? boundaryvalue(simplex, coc) : missing
            end
            s = findall(t -> !ismissing(t) && !iszero(t), m)
            if isempty(s)
                row.npts > length(cocycles) && continue
                model = cplx[row.npts]
                push!(cocycles[row.npts], (cocycle = onecochain(T, model, row.index), span = Span{U}(row.weight, U(Inf))))
            else
                li = last(s)
                lc, ic = c_small[li]
                for ii in s[1:end-1]
                    ratio = m[ii] / m[li]
                    c_small[ii].cocycle.values .-= ratio .* lc.values
                end
                (t0, t1) = endpoints(ic)
                c_small.span[li] = Span{U}(t0, row.weight)
            end
        end
    end
    map(cocycles) do c
        filter(t -> !isempty(t.span), c)
    end
end

persistent_cocycles(::Type{T}, cplx, sv, max_dim) where {T} =
    persistent_cocycles(T, eltype(cplx[1].values), cplx, sv, max_dim)

