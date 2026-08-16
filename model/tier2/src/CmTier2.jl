"""
    CmTier2

Tier-2 mechanistic cardiomyocyte cell-cycle fate model.

Two tiers over one shared node vocabulary (see `../MODEL.md` §"Model architecture"):

  * **Tier 1** (`../cmcycle`, Python) — a 55-node normalized-Hill logic network.
    Calibrated, but `tau` there is a relaxation constant, not a duration, so it
    cannot express absolute phase durations, duration *distributions*, ploidy and
    cell counting, cumulative-EdU versus instantaneous Ki-67, or FUCCI trace shape.

  * **Tier 2** (this package) — a mass-action / Michaelis-Menten ODE with real
    time in it, built on the published generic cell-cycle model in
    `Cell_Cycle_Model` (Womack, MIT). Those five questions are what it exists for.

# Phase 0 contract

This module currently exposes **only the inherited core, unchanged**. `src/inherited/`
is a byte-identical copy of the published `model_files/`; the golden regression test in
`test/` pins its behaviour. Cardiomyocyte modules (the E2F sub-family split, the
Ect2/RhoA/midbody cytokinesis arm, the maturation axis, nuclei-and-ploidy bookkeeping)
land in later phases and must reduce *exactly* to this baseline when disabled.

Do not edit `src/inherited/*.jl`. Extensions belong in their own files so the diff
against the published model stays legible to a reviewer.
"""
module CmTier2

using ComponentArrays
using DifferentialEquations
using Statistics
using FindPeaks1D

export state, state_names, params, modelDiffEq!
export solve_baseline, solve_drug, doubling_time, peak_period, species_index
export PUBLISHED_ALPHA, PUBLISHED_DOUBLING_TIME, PARAMFILE_ALPHA
# observables
export total_pools, fucci_state, fucci_fractions, phase_times
export FUCCI_THRESHOLD, PUBLISHED_FUCCI_THRESHOLD
# events
export EventThresholds, EventLog, solve_with_events, landmark_callbacks, trim
# fates
export Cycle, Bookkeeping, classify_cycles, quiescent, bookkeep, fate_summary, FATES_PHASE1
# tier 2 model
export tier2_state, tier2_state_names, tier2_params, tier2DiffEq!, solve_tier2
export e2f_repression, TIER2_SPECIES, TIER2_ENABLE_PARAMS

# ---------------------------------------------------------------------------
# The inherited core, verbatim.
#
# `include_all_scripts()` in the source repo glob-includes every .jl in a folder,
# which is how a stale second `modelDiffEq!` could silently win. Here the three
# files are named explicitly and `model_files/diff_eqn_plk1p_test.jl` (a 442-line
# shadow copy of diff_eqns.jl with different Km wiring) is deliberately not carried
# over.
# ---------------------------------------------------------------------------
include("inherited/parameters.jl")
include("inherited/state.jl")
include("inherited/diff_eqns.jl")

# ---------------------------------------------------------------------------
# The `α` provenance problem.
#
# The source repo carries FOUR different global time-scales, and the one in
# `parameters.jl` is not the one the published figures were made with. Because α
# multiplies the entire RHS, the cycle period scales exactly as 1/α, so this is a pure
# rescaling of every timing result:
#
#   α = 1.447  -> 28.11 h   hardcoded in `publication_fig_12_*_test.jl` and
#                           `*_varyingABE.jl`; reproduces the committed
#                           `experimental_vs_simulation_summary.csv` (DMSO 28.1 h)
#   α = 1.60   -> 25.43 h   `publication_fig_12_*_varyingVolo.jl`
#   α = 1.725  -> 23.58 h   implied by the default target in
#                           `Helpers/estimate_alpha_for_target_dt.jl`
#   α = 2.3    -> 17.69 h   the default in `parameters.jl` — matches no published number
#
# Tier 2 takes DEFAULT_ALPHA = the value the published figures actually used, and pins
# all four in the golden test so the choice stays visible. Verified: the pre- and
# post-`9ab1b1d` model files give identical periods, so the Km rewiring in that commit
# is numerically inert and is not the source of the spread.
# ---------------------------------------------------------------------------

"""Global time-scale used by the published figure scripts (`params(α = 1.447)`)."""
const PUBLISHED_ALPHA = 1.447

"""Doubling time the published figures report for DMSO control."""
const PUBLISHED_DOUBLING_TIME = 28.1

"""The `α` sitting in `parameters.jl`. Reproduces no published figure — see above."""
const PARAMFILE_ALPHA = 2.3

"""
    species_index(name) -> Int

Index of a species in the state vector. Throws rather than returning `nothing`, so a
typo fails loudly instead of silently indexing the wrong row.

Also resolves the Tier-2 species appended after the inherited 63, so observables, events
and fate code work unchanged against both `solve_baseline` and `solve_tier2` solutions —
the leading 63 indices are identical by construction.
"""
function species_index(name::AbstractString)
    idx = findfirst(==(name), state_names())
    idx === nothing || return idx
    extra = findfirst(==(name), TIER2_SPECIES)
    extra === nothing && error("unknown species $(repr(name)); see tier2_state_names()")
    return length(state_names()) + extra
end

"""
    solve_baseline(; alpha, tspan, con_VOL, con_ABE, kwargs...)

Solve the inherited model with no drug, using the same integrator the published
scripts use (`AutoTsit5(Rosenbrock23())`).

Returns the `ODESolution`. Tolerances are passed through so the golden test can pin
them; the published scripts rely on solver defaults, so the defaults here are those.
"""
function solve_baseline(; alpha::Real = PUBLISHED_ALPHA,
                          tspan::Tuple{<:Real,<:Real} = (0.0, 2500.0),
                          con_VOL::Real = 0.0,
                          con_ABE::Real = 0.0,
                          kwargs...)
    u0 = state()
    p = copy(params())
    p.α = alpha
    p.con_VOL = con_VOL
    p.con_ABE = con_ABE
    prob = ODEProblem(modelDiffEq!, u0, (Float64(tspan[1]), Float64(tspan[2])), p)
    return solve(prob, AutoTsit5(Rosenbrock23()); kwargs...)
end

"""
    doubling_time(sol; threshold, window) -> Float64

Mean cell-cycle period, in hours, from rising-edge crossings of `CCNB_CDK1` at
`threshold` × its windowed maximum, linearly interpolated between samples.

This is a faithful port of `doubling_time_for_alpha` in the source repo's
`Helpers/estimate_alpha_for_target_dt.jl` — the function that defined `α = 2.3`.
The measurement must not drift, or the golden test stops meaning anything.

Returns `NaN` when fewer than two crossings are found, matching the original.
"""
function doubling_time(sol; threshold::Float64 = 0.5,
                            window::Tuple{<:Real,<:Real} = (1800.0, 2400.0))
    idx = species_index("CCNB_CDK1")
    t = sol.t
    y = sol[idx, :]

    sel = findall(ti -> window[1] <= ti <= window[2], t)
    length(sel) < 3 && return NaN
    t_sub, y_sub = t[sel], y[sel]

    y_max = maximum(y_sub)
    y_max > 0 || return NaN
    y_norm = y_sub ./ y_max

    crossings = Float64[]
    for i in 2:length(y_norm)
        if y_norm[i-1] < threshold <= y_norm[i]
            t1, t2 = t_sub[i-1], t_sub[i]
            y1, y2 = y_norm[i-1], y_norm[i]
            push!(crossings, y2 != y1 ? t1 + (threshold - y1) / (y2 - y1) * (t2 - t1)
                                      : t_sub[i])
        end
    end

    length(crossings) < 2 && return NaN
    return mean(diff(crossings))
end

"""
    solve_drug(; vola_nM, abe_uM, alpha, tspan, t_add_drug)

Solve with volasertib (PLK1i, nM) and abemaciclib (CDK4/6i, µM) added at `t_add_drug`
via a `DiscreteCallback` that mutates `integrator.p` — the mechanism the published
`publication_fig_12_*` scripts use.

Note the callback condition is `t > t_add_drug`, i.e. it re-fires on every step after
the dose. That is idempotent here (it assigns constants) and is reproduced rather than
"fixed", because Phase 0's job is fidelity.
"""
function solve_drug(; vola_nM::Real = 0.0,
                      abe_uM::Real = 0.0,
                      alpha::Real = PUBLISHED_ALPHA,
                      tspan::Tuple{<:Real,<:Real} = (0.0, 2500.0),
                      t_add_drug::Real = 2000.0,
                      kwargs...)
    u0 = state()
    p = copy(params())
    p.α = alpha
    p.con_VOL = 0.0
    p.con_ABE = 0.0

    affect!(integrator) = begin
        integrator.p.con_VOL = Float64(vola_nM)
        integrator.p.con_ABE = Float64(abe_uM)
    end
    cb = DiscreteCallback((u, t, integrator) -> t > t_add_drug, affect!)

    prob = ODEProblem(modelDiffEq!, u0, (Float64(tspan[1]), Float64(tspan[2])), p)
    return solve(prob, AutoTsit5(Rosenbrock23()); callback = cb, kwargs...)
end

"""
    peak_period(sol, species="CCNB_CDK1"; window, height, prominence) -> Float64

Mean peak-to-peak interval, the measurement the published figure scripts use
(`find_peak_stats` with `height = prominence = 0.2`, then `mean(diff(peak_times))`).

Kept alongside [`doubling_time`](@ref), which uses the threshold-crossing method from
`Helpers/estimate_alpha_for_target_dt.jl`. The two are independent implementations of
the same quantity and agree to ~2 ms on the baseline oscillation; the golden test
asserts that, so a regression in either is caught.
"""
function peak_period(sol, species::AbstractString = "CCNB_CDK1";
                     window::Tuple{<:Real,<:Real} = (1800.0, 2000.0),
                     height::Real = 0.2, prominence::Real = 0.2)
    idx = species_index(species)
    t, y = sol.t, sol[idx, :]
    sel = findall(ti -> window[1] <= ti <= window[2], t)
    length(sel) < 2 && return NaN
    t_sub, y_sub = t[sel], y[sel]
    peaks, _ = findpeaks1d(y_sub; height = height, prominence = prominence)
    length(peaks) < 2 && return NaN
    return mean(diff(t_sub[peaks]))
end

# ---------------------------------------------------------------------------
# Phase 1: observables, discrete landmarks, fate classification.
# These add no biology and no parameters — they only measure the inherited model.
# ---------------------------------------------------------------------------
include("observables.jl")
include("events.jl")
include("fates.jl")

# ---------------------------------------------------------------------------
# Phase 2: the cardiomyocyte modules. Additive corrections over the inherited RHS,
# each identically zero at default parameters — see tier2_model.jl.
# ---------------------------------------------------------------------------
include("tier2_model.jl")

end # module
