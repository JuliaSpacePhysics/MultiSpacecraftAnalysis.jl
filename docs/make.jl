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
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "MMS comparison" => "mms_comparison.md",
    ],
    plugins = [bib],

)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/MultiSpacecraftAnalysis.jl",
    push_preview = true,
)
