using Documenter, MetabolicEP

push!(LOAD_PATH, "../src/")
makedocs(
    sitename = "MetabolicEP",
    modules = [MetabolicEP],
    doctest = false)
deploydocs(
    branch = "gh-pages",
    repo = "github.com/anna-pa-m/Metabolic-EP",
    versions = ["stable" => "v^"],
)
