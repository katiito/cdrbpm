#!/usr/bin/env julia
# =============================================================================
# compute_pairwise_distances.jl
# =============================================================================
# Computes pairwise genetic distances between all sequenced individuals within
# each simulation, mimicking what HIV-TRACE would produce from aligned sequences.
#
# Background:
#   The main simulation (j0.jl) only stores genetic distances between direct
#   donor-recipient pairs. HIV-TRACE computes distances between ALL sequenced
#   individuals, including those not directly linked by transmission (e.g.,
#   "siblings" both infected by the same source). This script reconstructs
#   those pairwise distances from the transmission tree.
#
# Method:
#   For any pair of sequenced individuals (u, v) with most recent common
#   ancestor (MRCA) individual m:
#     - Find the child of m on the path to u (cu) and to v (cv)
#     - t_mrca = min(cu.timeinfected, cv.timeinfected)  [the split point]
#     - total_evo_time = (u.timesequenced - t_mrca) + (v.timesequenced - t_mrca)
#     - distance ~ NegBin(s*mu*total_evo_time/omega, 1/(1+omega)) / s
#   This is consistent with the direct-pair formula in j0.jl and standard
#   molecular clock phylogenetics.
#
# Usage (CLI):
#   julia src/compute_pairwise_distances.jl <G_file> <D_file> [threshold] [out_prefix]
#
#   julia src/compute_pairwise_distances.jl \
#       src/experiment1-N50000-gens7-G.csv \
#       src/experiment1-N50000-gens7-D.csv \
#       0.015 \
#       src/experiment1-N50000-gens7
#
# Outputs:
#   <out_prefix>-pairs.csv   — all pairs with distance <= threshold
#     columns: simid, id1, id2, distance
#   <out_prefix>-clusters.csv — cluster membership for each sequenced individual
#     columns: simid, pid, cluster_id, cluster_size
#
# Parameters (matching j0.jl defaults):
#   mu    = 0.001/365   substitutions/site/day (mean clock rate)
#   omega = 0.5         variance inflation (overdispersion)
#   s     = 1000        precision scalar (sites)
# =============================================================================

try
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
catch err
    @warn "Unable to activate project" exception=(err, catch_backtrace())
end

using CSV
using DataFrames
using Distributions
using Random
using Printf

# -----------------------------------------------------------------------------
# Parameters (must match j0.jl)
# -----------------------------------------------------------------------------
const MU    = 0.001 / 365.0   # substitutions/site/day
const OMEGA = 0.5              # overdispersion
const S     = 1000             # precision scalar

# -----------------------------------------------------------------------------
# Genetic distance between two individuals given total evolutionary time
# Uses the same NegBin model as simgendist() in j0.jl
# -----------------------------------------------------------------------------
function sample_gendist(total_evo_time::Float64)::Float64
    total_evo_time <= 0.0 && return 0.0
    r = S * MU * total_evo_time / OMEGA
    p = 1.0 / (1.0 + OMEGA)
    muts = rand(NegativeBinomial(r, p))
    return Float64(muts) / S
end

# -----------------------------------------------------------------------------
# Parse pid string into a vector of integer segments
# "0.1.2" -> [0, 1, 2]
# -----------------------------------------------------------------------------
function parse_pid(pid::AbstractString)::Vector{Int}
    parse.(Int, split(pid, '.'))
end

# -----------------------------------------------------------------------------
# Find MRCA of two pids and return:
#   (mrca_pid_str, child_of_mrca_towards_u, child_of_mrca_towards_v)
# If one is the ancestor of the other, child_towards_ancestor = nothing
# -----------------------------------------------------------------------------
function find_mrca(pid_u::Vector{Int}, pid_v::Vector{Int})
    n = min(length(pid_u), length(pid_v))
    common_len = 0
    for i in 1:n
        pid_u[i] == pid_v[i] || break
        common_len = i
    end

    mrca_parts = pid_u[1:common_len]
    mrca_str   = join(mrca_parts, '.')

    # Child of MRCA on path to u (nothing if u IS the MRCA)
    cu_str = common_len < length(pid_u) ?
        join(pid_u[1:common_len+1], '.') : nothing

    # Child of MRCA on path to v (nothing if v IS the MRCA)
    cv_str = common_len < length(pid_v) ?
        join(pid_v[1:common_len+1], '.') : nothing

    (mrca_str, cu_str, cv_str)
end

# -----------------------------------------------------------------------------
# Build a lookup: pid_str -> timeinfected  (for finding t_mrca)
# -----------------------------------------------------------------------------
function build_tinf_lookup(G::DataFrame)::Dict{String,Float64}
    Dict(row.pid => row.timeinfected for row in eachrow(G))
end

# -----------------------------------------------------------------------------
# Compute total evolutionary time for a pair (u, v)
# Returns Inf if MRCA timeinfected is not found (should not happen in a
# well-formed tree, but can occur if ancestral individuals are unsequenced
# and their pid is not in G — we look up in the transmission-derived lookup)
# -----------------------------------------------------------------------------
function total_evo_time(
    pid_u::String, tseq_u::Float64,
    pid_v::String, tseq_v::Float64,
    tinf_lookup::Dict{String,Float64}
)::Float64

    pu = parse_pid(pid_u)
    pv = parse_pid(pid_v)

    _, cu_str, cv_str = find_mrca(pu, pv)

    # Determine t_mrca: the transmission time that splits the two lineages.
    # This is min(timeinfected of cu, timeinfected of cv).
    # If one individual IS the MRCA (cu_str or cv_str is nothing), use the
    # other child's timeinfected as the split point.
    t_mrca = if isnothing(cu_str) && isnothing(cv_str)
        # u == v (shouldn't be called for identical pids)
        return 0.0
    elseif isnothing(cu_str)
        # u is the MRCA; split is when u transmitted to cv
        get(tinf_lookup, cv_str, Inf)
    elseif isnothing(cv_str)
        # v is the MRCA; split is when v transmitted to cu
        get(tinf_lookup, cu_str, Inf)
    else
        t_cu = get(tinf_lookup, cu_str, Inf)
        t_cv = get(tinf_lookup, cv_str, Inf)
        min(t_cu, t_cv)
    end

    isinf(t_mrca) && return Inf

    # Both branches evolve independently from t_mrca
    (tseq_u - t_mrca) + (tseq_v - t_mrca)
end

# -----------------------------------------------------------------------------
# Process one simulation: compute all pairwise distances among sequenced
# individuals and return edges below threshold
# -----------------------------------------------------------------------------
function process_sim(
    G::DataFrame,
    threshold::Float64,
    tinf_lookup::Dict{String,Float64}
)
    # Only consider individuals who were sequenced
    Gseq = G[isfinite.(G.timesequenced), :]
    n = nrow(Gseq)

    pairs = Tuple{String,String,Float64}[]

    n < 2 && return pairs

    pids  = Gseq.pid
    tseqs = Gseq.timesequenced

    for i in 1:n-1
        for j in i+1:n
            t = total_evo_time(pids[i], tseqs[i], pids[j], tseqs[j], tinf_lookup)
            isinf(t) && continue
            t < 0.0  && continue    # can happen if sequenced before MRCA (numerical edge)
            d = sample_gendist(t)
            d <= threshold && push!(pairs, (pids[i], pids[j], d))
        end
    end

    pairs
end

# -----------------------------------------------------------------------------
# Find connected components (clusters) from edge list
# Returns: Dict pid -> cluster_id
# -----------------------------------------------------------------------------
function connected_components(
    all_pids::Vector{String},
    edges::Vector{Tuple{String,String,Float64}}
)::Dict{String,Int}

    # Union-find
    parent = Dict(p => p for p in all_pids)

    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]   # path compression
            x = parent[x]
        end
        x
    end

    function union!(a, b)
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end

    for (u, v, _) in edges
        haskey(parent, u) && haskey(parent, v) && union!(u, v)
    end

    # Assign integer cluster IDs
    roots  = Dict{String,Int}()
    result = Dict{String,Int}()
    cid    = 0
    for p in all_pids
        r = find(p)
        if !haskey(roots, r)
            cid += 1
            roots[r] = cid
        end
        result[p] = roots[r]
    end

    result
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function compute_pairwise_distances(
    g_file::String;
    threshold::Float64 = 0.015,
    out_prefix::String = "",
    progress_every::Int = 1000
)
    println("Loading data...")
    Gall = CSV.read(g_file, DataFrame; stringtype=String)
    Gall.simid = string.(Gall.simid)   # coerce to String in case CSV inferred numeric

    simids = unique(Gall.simid)
    n_sims = length(simids)
    println("Processing $n_sims simulations (threshold=$threshold)...")

    # Index rows by simid (avoids groupby ambiguity with Base.Iterators.groupby)
    sim_rows = Dict{String, Vector{Int}}()
    for i in 1:nrow(Gall)
        push!(get!(sim_rows, Gall.simid[i], Int[]), i)
    end

    t_start = time()

    # Accumulators
    all_pairs    = DataFrame(simid=String[], id1=String[], id2=String[], distance=Float64[])
    all_clusters = DataFrame(simid=String[], pid=String[], cluster_id=Int[], cluster_size=Int[])

    for (idx, simid) in enumerate(simids)
        G1 = Gall[sim_rows[simid], :]

        # Build tinf lookup restricted to this sim's pids
        tinf_lookup = Dict(row.pid => row.timeinfected for row in eachrow(G1))

        pairs = process_sim(G1, threshold, tinf_lookup)

        # Cluster membership for all sequenced individuals in this sim
        Gseq = G1[isfinite.(G1.timesequenced), :]
        if nrow(Gseq) > 0
            clust = connected_components(Gseq.pid, pairs)
            # Count cluster sizes
            sizes = Dict{Int,Int}()
            for cid in values(clust)
                sizes[cid] = get(sizes, cid, 0) + 1
            end
            for pid in Gseq.pid
                cid  = clust[pid]
                push!(all_clusters, (simid, pid, cid, sizes[cid]))
            end
        end

        # Append pairs
        for (u, v, d) in pairs
            push!(all_pairs, (simid, u, v, d))
        end

        if progress_every > 0 && (idx % progress_every == 0 || idx == n_sims)
            elapsed = time() - t_start
            rate    = idx / elapsed
            eta     = (n_sims - idx) / rate
            @printf("  [%6d / %d] %.1f sims/sec, elapsed %.1fs, ETA %.1fs\n",
                    idx, n_sims, rate, elapsed, eta)
        end
    end

    total_time = time() - t_start
    println("\nDone in $(round(total_time; digits=1))s")
    println("  Pairs below threshold: $(nrow(all_pairs))")
    println("  Sequenced individuals with cluster assignment: $(nrow(all_clusters))")

    # Write outputs
    if !isempty(out_prefix)
        pairs_path    = string(out_prefix, "-pairs.csv")
        clusters_path = string(out_prefix, "-clusters.csv")
        CSV.write(pairs_path, all_pairs)
        CSV.write(clusters_path, all_clusters)
        println("  Wrote: $pairs_path")
        println("  Wrote: $clusters_path")
    end

    (; pairs=all_pairs, clusters=all_clusters)
end

# -----------------------------------------------------------------------------
# CLI entrypoint
# -----------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    g_file     = length(ARGS) >= 1 ? ARGS[1] : "src/experiment1-N50000-gens7-G.csv"
    threshold  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.015
    out_prefix = length(ARGS) >= 3 ? ARGS[3] : replace(g_file, r"-G\.csv$" => "")

    compute_pairwise_distances(g_file; threshold=threshold, out_prefix=out_prefix)
end
