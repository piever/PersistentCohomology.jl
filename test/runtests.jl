using PersistentCohomology
using PersistentCohomology: Span

using Test

using SparseArrays, StructArrays
using LinearAlgebra
using GaloisFields

@testset "vietorisrips" begin

    m = [0.0 0.2 0.9 0.0
         0.2 0.0 0.5 0.0
         0.9 0.5 0.0 0.0
         0.0 0.0 0.0 0.0]

    s = sparse(m)

    cplx = vietorisrips(s, 2)
    @test keys(cplx[1]) == StructArray([(1,), (2,), (3,), (4,)],)
    @test values(cplx[1]) == [0.0, 0.0, 0.0, 0.0]

    @test keys(cplx[2]) == StructArray([(2, 1), (3, 1), (3, 2)])
    @test values(cplx[2]) == [0.2, 0.9, 0.5]

    @test keys(cplx[3]) == StructArray([(3, 2, 1)])
    @test values(cplx[3]) == [0.9]

    @test length(cplx) == 3

    cplx = vietorisrips(s, 1)
    @test keys(cplx[1]) == StructArray([(1,), (2,), (3,), (4,)],)
    @test values(cplx[1]) == [0.0, 0.0, 0.0, 0.0]

    @test keys(cplx[2]) == StructArray([(2, 1), (3, 1), (3, 2)])
    @test values(cplx[2]) == [0.2, 0.9, 0.5]

    @test length(cplx) == 2
end
    
@testset "persistent_cocycles" begin
    GF = @GaloisField 3
    ts = [-pi, -pi/2, 0, pi/2]
    pts = map(t -> [cos(t), sin(t)], ts)

    m = round.([norm(pt1-pt2) for pt1 in pts, pt2 in pts], digits = 4)
    cplx = vietorisrips(sparse(m), 2)
    @test keys(cplx[1]) == StructArray(([1, 2, 3, 4],))
    @test keys(cplx[2]) == StructArray([(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)])
    @test keys(cplx[3]) == StructArray([(3, 2, 1), (4, 2, 1), (4, 3, 1), (4, 3, 2)])

    short = m[1, 2]
    long = m[1, 3]

    cocycles, spans = persistent_cocycles(GF, cplx, 1)

    @test spans[1][1] == Span{Float64}(0, Inf)
    @test findnz(values(cocycles[1][1])) == ([1, 2, 3, 4], GF.([1, 1, 1, 1]))
    @test spans[1][2] == Span{Float64}(0, short)
    @test findnz(values(cocycles[1][2])) == ([2], GF.([1]))
    @test spans[1][3] == Span{Float64}(0, short)
    @test findnz(values(cocycles[1][3])) == ([3], GF.([1]))
    @test spans[1][4] == Span{Float64}(0, short)
    @test findnz(values(cocycles[1][4])) == ([4], GF.([1]))
    @test length(cocycles[1]) == 4

    @test spans[2][1] == Span{Float64}(short, long)
    @test findnz(values(cocycles[2][1])) == ([6], GF.([1]))
    @test length(cocycles[2]) == 1
end
