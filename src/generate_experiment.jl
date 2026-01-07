#!/usr/bin/env julia

# Batch generator: run multiple simulations and write combined CSVs
# Usage examples:
#   julia src/generate_experiment.jl 100 experiment1 5
#     -> runs 100 sims, maxgenerations=5, writes src/experiment1-N100-D.csv and -G.csv
#   julia src/generate_experiment.jl 50 myexp
#     -> runs 50 sims, maxgenerations=5 (default), writes src/myexp-N50-*.csv

include("j0.jl")
using CSV
using DataFrames

function generate_experiment(n_sims::Int;
    out_prefix::AbstractString = "experiment1",
    maxgenerations::Int = 5,
    initialcontact::Symbol = :G,
)
    @assert n_sims > 0 "n_sims must be positive"
    G_all = DataFrame()
    D_all = DataFrame()

    for i in 1:n_sims
        o = simbp(P; maxgenerations=maxgenerations, initialcontact=initialcontact)
        if i == 1
            G_all = copy(o.G)
            D_all = copy(o.D)
        else
            G_all = vcat(G_all, o.G; cols=:union)
            D_all = vcat(D_all, o.D; cols=:union)
        end
    end

    # Build output filenames with N tag
    tag = string(out_prefix, "-N", n_sims)
    out_dir = @__DIR__  # write into src/ so the R script finds it via here('src', ...)
    g_path = joinpath(out_dir, string(tag, "-G.csv"))
    d_path = joinpath(out_dir, string(tag, "-D.csv"))

    CSV.write(g_path, G_all)
    CSV.write(d_path, D_all)

    println("Wrote: \n  ", g_path, " (", nrow(G_all), " rows)\n  ", d_path, " (", nrow(D_all), " rows)")
    return (; G=G_all, D=D_all, g_path, d_path)
end

# CLI entrypoint
if abspath(PROGRAM_FILE) == @__FILE__
    n_sims = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    out_prefix = length(ARGS) >= 2 ? ARGS[2] : "experiment1"
    maxgenerations = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 5
    generate_experiment(n_sims; out_prefix=out_prefix, maxgenerations=maxgenerations)
end
