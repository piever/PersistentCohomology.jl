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

    @test cplx[2].simplices == StructArray([(2, 1), (3, 1), (3, 2)])
    @test cplx[2].values == [0.2, 0.9, 0.5]

    @test cplx[3].simplices == StructArray([(3, 2, 1)])
    @test cplx[3].values == [0.9]

    @test length(cplx) == 3

    cplx = vietorisrips(s, 1)
    @test cplx[1].simplices == StructArray([(1,), (2,), (3,), (4,)],)
    @test cplx[1].values == [0.0, 0.0, 0.0, 0.0]

    @test cplx[2].simplices == StructArray([(2, 1), (3, 1), (3, 2)])
    @test cplx[2].values == [0.2, 0.9, 0.5]

    @test length(cplx) == 2
end
    
# using GaloisFields
# const GF = @GaloisField 3
# ts = [-pi, -pi/2, 0, pi/2]
# pts = map(t -> [cos(t), sin(t)], ts)

# m = [norm(pt1-pt2) for pt1 in pts, pt2 in pts]
# cplx = vietorisrips(sparse(m), 2)
# sv = PersistentCohomology.simplexlist(cplx)

# PersistentCohomology.persistent_cocycles(GF, cplx, sv, 1)

# s = spzeros(4, 4)
# s[1, 2] = 0.3
# s[2, 3] = 0.4
# s[3, 4] = 0.5
# s[1, 4] = 0.3
# s[1, 3] = 0.7
# s[2, 4] = 0.8



# cplx = vietorisrips(s, 2)
# sv = PersistentCohomology.simplexlist(cplx)

# PersistentCohomology.persistent_cocycles(F, cplx, sv, 1)

