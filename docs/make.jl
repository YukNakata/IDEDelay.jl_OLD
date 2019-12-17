using Documenter, IDEDelay

makedocs(;
    modules=[IDEDelay],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/YukNakata/IDEDelay.jl/blob/{commit}{path}#L{line}",
    sitename="IDEDelay.jl",
    authors="YukNakata <yunayuna.na@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/YukNakata/IDEDelay.jl",
)
