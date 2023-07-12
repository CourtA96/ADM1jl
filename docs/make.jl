# push!(LOAD_PATH,"C:\\Users\\court\\OneDrive - University of Guelph\\PhD\\Research\\Code\\ADM1jl\\ADM1jl\\src")
#test


using Documenter, ADM1

makedocs(
    sitename="ADM1.jl",
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