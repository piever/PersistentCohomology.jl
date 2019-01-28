module PersistentCohomology

using StructArrays
using StructArrays: fieldarrays
using SparseArrays, LinearAlgebra
using IntervalSets

export Cochain, vietorisrips

include("cochain.jl")
include("vietorisrips.jl")
include("cocycle.jl")

end # module
