function emptysimplexvector(n)
    typ = NTuple{n+1, Int}
    StructVector{typ}(undef, 0)
end

function lower_neighbors(g, u)
    Iterators.filter(i -> i < u, neighbors(g, u))
end

function vietorisrips(g::AbstractSimpleWeightedGraph, max_dim = 2)
    v = Dict{Int, Cochain}()
    us = vertices(g)
    # v[0] = Cochain(StructVector((vs,)), fill(0.0, length(us)))
    for i in 0:max_dim
        v[i] = Cochain(emptysimplexvector(i), Float64[])
    end
    for u in us
        ln = lower_neighbors(g, u)
        add_cofaces!(g, max_dim, (u,), ln, v)
    end
    return v
end

function add_cofaces!(g, max_dim, sim, ln, v)
    push!(v[length(sim)-1].simplices, sim)
    if length(sim) > max_dim
        return
    else
        for n in ln
            ln1 = lower_neighbors(g, n)
            mn = Iterators.filter(i -> (i in ln1), ln)
            add_cofaces!(g, max_dim, (n, sim...) , mn, v)
        end
    end
end
