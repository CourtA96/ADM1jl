using Pkg

Pkg.add("Documenter")

using Documenter, ADM1jl

makedocs(
    sitename="ADM1jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "page2.md",
            "page3.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/CourtA96/ADM1jl.git",
    push_preview = true,
)