function initialize_cocycles(::Type{T}, ::Type{U}, cplx) where {T, U}
    instances = map(t -> Cochain(t, spzeros(T, length(t))), cplx)
    map(instances) do inst
        typ = Tuple{typeof(inst), Interval{:closed, :open, U}}
        StructVector{typ}(undef, 0)
    end
end

function persistent_cocycles(::Type{T}, ::Type{U}, cplx, sv, max_dim) where {T, U}
    cocycles = initialize_cocycles(T, U, cplx)
    for row in sv
        row.npts == 1 && continue
        c_small = cocycles[row.npts-1]
        simplex = cplx[row.npts].simplices[row.index]
        m = map(c_small) do (coc, interval)
            row.weight in interval ? boundaryvalue(simplex, coc) : missing
        end
        s = findall(t -> !ismissing(t) && !iszero(t), m)
        if isempty(s)
            row.npts <= max_dim && continue
            model = cplx[row.npts]
            null = spzeros(T, length(model))
            null[row.index] = one(T)
            push!(cocycles[row.npts], (Cochain(model, null), Interval{:closed, :open, U}(row.weight, U(Inf))))
        else
            li = last(s)
            lc, ic = c_small[li]
            for ii in s[1:end-1]
                ratio = m[ii] / m[li]
                c_small[ii][1].values .-= ratio .* lc.values
            end
            (t0, t1) = endpoints(ic)
            fieldarrays(c_small)[2][li] = Interval{:closed, :open, U}(t0, row.weight)
        end
    end
    cocycles
end
