module PersistentCohomology

using StructArrays
using StructArrays: finduniquesorted
using SparseArrays, LinearAlgebra
using IntervalSets

export Cochain, vietorisrips, persistent_cocycles

include("cochain.jl")
include("vietorisrips.jl")
include("cocycle.jl")
include("itertools.jl")

end # module
