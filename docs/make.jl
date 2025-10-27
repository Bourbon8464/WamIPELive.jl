# docs/make.jl
using Documenter
using WamIPELive

cd(@__DIR__)

makedocs(
    modules   = [WamIPELive],
    sitename  = "SpaceAGORA.jl",
    format    = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
    pagesonly = true,
    checkdocs = :exports,             # or :none if you truly don’t want checks
    warnonly  = [:missing_docs, :cross_references],
    pages = [
        "Introduction" => "index.md",
        "Packages" => Any[
            "Custom Packages" => "packages/introduction.md",
            "WamIPELive.jl" => "packages/WamIPELive.md",
            "WamIPEDensity.jl" => "packages/WamIPEDensity.md",
        ],
    ],
)

deploydocs(; repo = "github.com/Bourbon8464/WamIPELive.jl",  
branch = "gh-pages", devbranch = "main")
