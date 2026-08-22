# Derived observables — the quantities the two FUCCI papers actually measure.
#
# Tier 1 cannot produce any of these (see ../MODEL.md: "What Tier 1 deliberately cannot
# do"). They are the reason Tier 2 exists, and they are what it will be calibrated
# against, using the target rows Tier 1 structurally could not use.

"""
    total_pools(sol) -> NamedTuple

The summed pools the published `post_process_cc` adds to a solution: `WEE1T`, `PLK1T`,
`APCCT`, and `GemininT`.

`GemininT` matters and was **commented out** in the source repo's `post_process_cc`,
while `plot_fucci_backgrounds` reads `solution_df.GemininT` — so the FUCCI shading path
there would raise. It is restored here because total geminin is the green FUCCI channel:
the reporter reports geminin whether or not it is bound to Cdt1.
"""
function total_pools(sol)
    g(s) = sol[species_index(s), :]
    return (WEE1T   = g("WEE1")    .+ g("WEE1p"),
            PLK1T   = g("PLK1")    .+ g("PLK1p"),
            APCCT   = g("APCC")    .+ g("APCCP"),
            GemininT = g("Geminin") .+ g("Geminin_CDT1"))
end

"""
Threshold hardcoded in the source repo's `plot_fucci_backgrounds`. It is a **plotting
default**, not a calibrated value — see [`FUCCI_THRESHOLD`](@ref).
"""
const PUBLISHED_FUCCI_THRESHOLD = 0.05

"""
    FUCCI_THRESHOLD

Reporter positivity cutoff, **calibrated: 0.02**.

This is a property of the *reporter and the microscope*, not of the model — it is where a
fluorescence intensity is called positive — so leaving it at an arbitrary value and then
concluding things about the biology is a mistake. It is counted as one declared fitted
parameter.

## Why 0.02, and how it was set without touching any published measurement

The model reports phase durations two independent ways: from cyclin peaks
([`phase_times`](@ref), which never looks at Cdt1 or geminin) and from the FUCCI
channels ([`fucci_fractions`](@ref), which never looks at a cyclin). Requiring the two to
agree fixes the cutoff with no external data at all. Scanning it:

| cutoff | Cdt1-positive | cyclin G1 | geminin-positive | cyclin S+G2+M | double-negative |
|---|---|---|---|---|---|
| 0.050 | 0.517 | 0.613 | 0.127 | 0.387 | **0.356** |
| 0.025 | 0.581 | 0.613 | 0.314 | 0.387 | 0.105 |
| **0.020** | **0.609** | **0.613** | **0.391** | **0.387** | **0.000** |
| 0.015 | 0.813 | 0.617 | 0.187 | 0.387 | 0.000 |

A sharp optimum: total error 0.0006 at 0.020 against 0.211 at 0.025 and 0.400 at 0.015.
Two measurements that share no equations agree to under one percentage point, and the
double-negative population — which a real asynchronous FUCCI culture does not have —
vanishes exactly there.

**This retracts the Phase 1 "the FUCCI layer has no G1/S state" finding.** That was an
artefact of the 0.05 plotting default, not a defect: at 0.05, geminin clears the cutoff
for only 12.6 % of the cycle and never overlaps Cdt1, so the double-positive is empty and
36 % of the cycle reads double-negative. The licensing layer was fine; the ruler was
wrong. The inherited Cdt1/geminin dynamics reproduce the model's own phase structure,
which is a genuine internal validation and the reason this model is a better Tier-2 base
than one without FUCCI observables at all.
"""
const FUCCI_THRESHOLD = 0.02

"""
    fucci_state(cdt1, geminin_total; threshold=FUCCI_THRESHOLD) -> Symbol

The four-way FUCCI classification both papers report.

| Cdt1 (red) | geminin (green) | state | reported as |
|---|---|---|---|
| + | − | `:G0G1` | mKO2-only |
| + | + | `:G1S`  | double-positive |
| − | + | `:SG2M` | mAG-only |
| − | − | `:negative` | double-negative |

`:G1S` being a *double-positive* rather than a transition instant is what makes
Baniol's Suppl 1G correction possible — only 22.7 % of P0 double-positives are Ki-67
positive, the rest having prematurely exited. That correction is what turned Tier 1's
mouse validation into a 1.03x match (`../TODO.md` item 6).
"""
function fucci_state(cdt1::Real, geminin_total::Real; threshold::Real = FUCCI_THRESHOLD)
    red   = cdt1 > threshold
    green = geminin_total > threshold
    red && green  && return :G1S
    red           && return :G0G1
    green         && return :SG2M
    return :negative
end

"""
    fucci_fractions(sol; window, threshold) -> Dict{Symbol,Float64}

Time-weighted fraction of `window` spent in each FUCCI state.

For an asynchronous population at steady state, the fraction of *time* one cell spends
in a state equals the fraction of *cells* found in it — which is what flow cytometry
measures. That equivalence is what makes this directly comparable to the ten
`fucci_fraction` rows in `../cmcycle/data/cmcycle_targets.csv` (Murganti Fig 1C,
Baniol Fig 1D).

Trapezoidal in time, so it does not inherit the solver's non-uniform step density as a
weighting — a plain `count`/`length` over saved points would.
"""
function fucci_fractions(sol; window::Tuple{<:Real,<:Real} = (1800.0, 2500.0),
                              threshold::Real = FUCCI_THRESHOLD)
    pools = total_pools(sol)
    cdt1  = sol[species_index("CDT1"), :]
    gemT  = pools.GemininT
    t     = sol.t

    sel = findall(ti -> window[1] <= ti <= window[2], t)
    length(sel) < 2 && error("window $(window) contains < 2 saved points")

    acc = Dict(:G0G1 => 0.0, :G1S => 0.0, :SG2M => 0.0, :negative => 0.0)
    total = 0.0
    for k in 1:(length(sel) - 1)
        i, j = sel[k], sel[k+1]
        dt = t[j] - t[i]
        dt <= 0 && continue
        # midpoint state over the interval
        st = fucci_state((cdt1[i] + cdt1[j]) / 2, (gemT[i] + gemT[j]) / 2;
                         threshold = threshold)
        acc[st] += dt
        total   += dt
    end
    total > 0 || error("window $(window) has zero duration")
    return Dict(k => v / total for (k, v) in acc)
end

"""
    phase_times(sol; window, mito_end_threshold) -> NamedTuple

Cell-cycle phase durations, using the published boundary definitions from
`Functions Julia/determine_cc_phase_times.jl`:

  * **G1/S** — CCNE peak
  * **S/G2** — PLK1p rising through half its peak value
  * **G2/M** — CDC25Cp peak
  * **M/G1** — CCNB_CDK1 falling below `mito_end_threshold`

Reimplemented rather than reused so it returns `nothing` cleanly on failure. The
original assigns `s_time`/`g2_time`/`m_time` only inside its `if success` branch but
references them unconditionally in its return, so an unsuccessful call raises
`UndefVarError` instead of reporting failure.

Returns `nothing` when the boundaries are missing or out of order — which is itself the
check that the model is behaving like a cell cycle.

## Peak criteria are relative, and have to be

The original takes a per-variable `peak_params` dict because a single absolute `height`
cannot serve all four species: at α = 1.447 free `CCNE` peaks at 0.132, so the
`height = 0.2` the figure scripts use for `CCNB_CDK1` (peak 0.662) finds **zero** CCNE
peaks and the G1/S boundary is never located. Defaults here are therefore a fraction of
each species' own windowed maximum, which is scale-free across species and across the
maturation range Phase 2 introduces. `rel_height`/`rel_prominence` override it.

Note this is a different convention from [`peak_period`](@ref) and
[`doubling_time`](@ref), which keep the published absolute 0.2 because they reproduce
published numbers and are pinned by test.
"""
function phase_times(sol; window::Tuple{<:Real,<:Real} = (1800.0, 2000.0),
                          mito_end_threshold::Real = 0.1,
                          rel_height::Real = 0.2, rel_prominence::Real = 0.2)
    t = sol.t
    sel = findall(ti -> window[1] <= ti <= window[2], t)
    length(sel) < 10 && return nothing
    ts = t[sel]
    grab(s) = sol[species_index(s), :][sel]

    ccne, plk1p, cdc25cp = grab("CCNE"), grab("PLK1p"), grab("CDC25Cp")
    mito = grab("CCNB_CDK1")

    function pk(v)
        m = maximum(v)
        m > 0 || return Int[]
        first(findpeaks1d(v; height = rel_height * m, prominence = rel_prominence * m))
    end
    ccne_pk, plk1p_pk, cdc25cp_pk, mito_pk = pk(ccne), pk(plk1p), pk(cdc25cp), pk(mito)
    (isempty(ccne_pk) || isempty(plk1p_pk) || isempty(cdc25cp_pk) || isempty(mito_pk)) &&
        return nothing

    g1s = ts[ccne_pk]
    g2m = ts[cdc25cp_pk]

    # S/G2: last index at or below half the peak, searching backwards from the peak.
    sg2_idx = [findlast(x -> x <= plk1p[m] / 2, plk1p[1:m]) for m in plk1p_pk]
    any(isnothing, sg2_idx) && return nothing
    sg2 = ts[Int.(sg2_idx)]

    # M/G1: first drop below threshold after each mitotic peak.
    mg1 = Float64[]
    for m in mito_pk
        i = findfirst(j -> mito[j] < mito_end_threshold, (m+1):length(mito))
        i === nothing || push!(mg1, ts[m + i])
    end

    (length(g1s) >= 2 && length(sg2) >= 2 && length(g2m) >= 2 && length(mg1) >= 1) ||
        return nothing
    (g1s[1] < sg2[1] < g2m[1] < mg1[1]) || return nothing

    cell_cycle = g1s[2] - g1s[1]
    return (cell_cycle = cell_cycle,
            g1 = g1s[2] - mg1[1],
            s  = sg2[1] - g1s[1],
            g2 = g2m[1] - sg2[1],
            m  = mg1[1] - g2m[1],
            g1_percent = 100 * (g1s[2] - mg1[1]) / cell_cycle,
            boundaries = (g1s = g1s, sg2 = sg2, g2m = g2m, mg1 = mg1))
end
