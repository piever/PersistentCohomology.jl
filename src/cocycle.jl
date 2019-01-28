function initialize_cocycles(::Type{T}, cplx) where {T}
    instances = map(t -> Cochain(t, spzeros(T, length(t))), cplx)
    map(instances) do inst
        typ = Tuple{typeof(inst), Interval{:closed, :open, eltype(cplx[1].values)}}
        StructVector{typ}(undef, 0)
    end
end
