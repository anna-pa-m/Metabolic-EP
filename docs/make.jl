using Documenter, MetabolicEP

push!(LOAD_PATH,"../src/")
makedocs(sitename = "",
	 modules = [MetabolicEP],
	 doctest = true)
deploydocs(
	   branch = "gh-pages",
	   repo = "github.com/anna-pa-m/Metabolic-EP",
	   versions = ["stable" => "v^"]
	  )
