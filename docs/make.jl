# push!(LOAD_PATH,"C:\\Users\\court\\OneDrive - University of Guelph\\PhD\\Research\\Code\\ADM1jl\\ADM1jl\\src")

using Documenter, ADM1code

makedocs(
    sitename="ADM1code.jl",
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
)