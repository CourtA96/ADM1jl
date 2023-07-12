# push!(LOAD_PATH,"C:\\Users\\court\\OneDrive - University of Guelph\\PhD\\Research\\Code\\ADM1jl\\ADM1jl\\src")

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
    forcepush = false,
    deploy_config = auto_detect_deploy_system(),
    push_preview = false,
    repo_previews = repo,
    branch_previews = branch,
)