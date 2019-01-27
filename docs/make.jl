using Documenter, PersistentCohomology

makedocs(
    modules = [PersistentCohomology],
    format = :html,
    sitename = "PersistentCohomology.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/piever/PersistentCohomology.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing,
)
