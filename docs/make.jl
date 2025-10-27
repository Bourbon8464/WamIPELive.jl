# docs/make.jl
using Documenter
using WamIPELive

cd(@__DIR__)

makedocs(
    modules   = [WamIPELive],
    sitename  = "WamIPELive.jl",
    format    = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
    pagesonly = true,
    checkdocs = :exports,             # or :none if you truly don’t want checks
    warnonly  = [:missing_docs, :cross_references],
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(; repo = "github.com/Bourbon8464/WamIPELive.jl.git", devbranch = "main")
