using PersistentCohomology
using Test

using SparseArrays, StructArrays

@testset "vietorisrips" begin

    m = [0.0 0.2 0.9 0.0
         0.2 0.0 0.5 0.0
         0.9 0.5 0.0 0.0
         0.0 0.0 0.0 0.0]

    s = sparse(m)

    cplx = vietorisrips(s, 2)
    @test cplx[1].simplices == StructArray([(1,), (2,), (3,), (4,)],)
    @test cplx[1].values == [0.0, 0.0, 0.0, 0.0]

    @test cplx[2].simplices == StructArray([(1, 2), (1, 3), (2, 3)])
    @test cplx[2].values == [0.2, 0.9, 0.5]

    @test cplx[3].simplices == StructArray([(1, 2, 3)])
    @test cplx[3].values == [0.9]

    @test length(cplx) == 3

    cplx = vietorisrips(s, 1)
    @test cplx[1].simplices == StructArray([(1,), (2,), (3,), (4,)],)
    @test cplx[1].values == [0.0, 0.0, 0.0, 0.0]

    @test cplx[2].simplices == StructArray([(1, 2), (1, 3), (2, 3)])
    @test cplx[2].values == [0.2, 0.9, 0.5]

    @test length(cplx) == 2
end
    

