using MultiSpacecraftAnalysis
using Documenter
using DocumenterCitations

const bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"))

DocMeta.setdocmeta!(MultiSpacecraftAnalysis, :DocTestSetup, :(using MultiSpacecraftAnalysis); recursive = true)

makedocs(;
    modules = [MultiSpacecraftAnalysis],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "MultiSpacecraftAnalysis.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSpacePhysics.github.io/MultiSpacecraftAnalysis.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
    ],
    plugins = [bib],

)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/MultiSpacecraftAnalysis.jl",
    devbranch = "main",
)
