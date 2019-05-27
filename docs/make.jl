using Documenter, MetabolicEP

push!(LOAD_PATH,"../src/")
makedocs(sitename = "MetabolicEP",
         modules = [MetabolicEP],
         doctest = true)
deploydocs(
           branch = "gh-pages",
           repo = "github.com/anna-pa-m/MetabolicEP.git",
           versions = ["stable" => "v^"]
          )
