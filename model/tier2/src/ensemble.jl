# The population layer — `../TODO.md` item 4.
#
# ## Why a population at all
#
# A deterministic trajectory returns exactly one fate. Murganti report four fractions
# (90.35 / 5.09 / 3.16 / 1.40 %), which are a property of a heterogeneous population, and
# TODO.md item 4 is explicit that "fate fractions are tail statistics — reproducing 1.40 %
# needs the *shape* of the parameter distribution, not its mean, so a mean-field treatment
# cannot get there".
#
# This item is *blocked* in Tier 1 on performance (item 5: a steady state is ~405 ms and
# the ensemble needs 1e4-1e5 runs). It is close to free here: `EnsembleThreads` over the
# machine's cores.
#
# ## The noise scale is measured, not chosen
#
# Heterogeneity has exactly one free number — the width of the lognormal. It is NOT fitted
# to the fate fractions, which are the hold-out. It is set by requiring the simulated
# S/G2/M duration CV to match Baniol's measured 0.265 (P0 mouse, 15.1 +/- 4.0 SD), exactly
# as TODO.md item 4 specifies: "calibrate the noise scale to the measured duration CV
# (0.265 for mouse P0), not to a guessed value."

"""
Step cap per cell.

Sized from measurement, not guessed: a healthy cell integrates its span in ~7,600 steps,
so 100,000 is >13x headroom. The cap matters because heterogeneity draws occasional
pathological parameter sets and `Threads.@threads` makes wall-clock the SLOWEST cell — at
a 2,000,000 cap a single bad draw stalled a 64-cell run past 13 minutes. A capped cell is
scored `:Failed` and reported in the return value, never silently dropped.
"""
const MAXITERS_PER_CELL = 100_000

"""Measured S/G2/M duration CV, Baniol P0 mouse: 4.0 / 15.1."""
const MEASURED_DURATION_CV = 4.0 / 15.1

"""
    heterogeneous_params(p0, sigma, rng) -> ComponentVector

One cell's parameters: `p0` with every synthesis rate multiplied by an independent
lognormal factor of log-scale `sigma`, median 1.

Synthesis rates only. Cell-to-cell variability in a clonal population is dominated by
differences in expression level, not by differences in the chemistry — a rate constant is
a property of the molecule, an abundance is a property of the cell. Perturbing binding
constants would be asserting that two cells' proteins behave differently, which is a much
stronger claim than the data supports.

`exp(sigma*z)` rather than `1 + sigma*z`: concentrations are non-negative and
multiplicatively distributed, and a normal factor would put weight on negative synthesis.
"""
function heterogeneous_params(p0, sigma::Real, rng)
    p = copy(p0)
    sigma <= 0 && return p
    for k in keys(p)
        ks = String(k)
        startswith(ks, "ks_") || continue
        p[k] = p[k] * exp(sigma * randn(rng))
    end
    return p
end

"""
    run_ensemble(; n, sigma, enable, alpha, tspan, window, seed, batch) -> NamedTuple

Simulate `n` heterogeneous cells and return their fate fractions.

Each cell is scored by the fate of its **first completed cycle** in `window`, or
`:Quiescent` if it completes none — which is how Murganti score a live-imaging movie: one
cell, one observed outcome. Taking a modal fate over many cycles would be scoring
something the experiment cannot see.

Runs on `EnsembleThreads`; start Julia with `-t auto` to use the machine.

Returns fractions, counts, the Monte-Carlo standard error on each fraction, and the
simulated duration CV so the noise calibration can be checked in the same call.
"""
function run_ensemble(; n::Int = 1000,
                        sigma::Real = 0.10,
                        enable::NamedTuple = NamedTuple(),
                        alpha::Real = PUBLISHED_ALPHA,
                        tspan::Tuple{<:Real,<:Real} = (0.0, 1600.0),
                        window::Tuple{<:Real,<:Real} = (1000.0, 1500.0),
                        thresholds::EventThresholds = EventThresholds(),
                        seed::Int = 20260816,
                        cytokinesis::Bool = true)
    p0 = tier2_params(; enable...)
    p0.α = alpha
    u0 = tier2_state()

    labels = cytokinesis ? (:Quiescent, :Division, :Binucleation, :Polyploidization) :
                           (:Quiescent, :MitoticCompletion, :Polyploidization)
    fates = Vector{Symbol}(undef, n)
    durations = Vector{Float64}(undef, n)
    fill!(durations, NaN)

    Threads.@threads for i in 1:n
        # Per-cell RNG seeded from the run seed and the index, so a result is
        # reproducible and independent of thread scheduling.
        rng = Random.Xoshiro(hash((seed, i)))
        p = heterogeneous_params(p0, sigma, rng)
        log = EventLog()
        prob = ODEProblem(tier2DiffEq!, u0, (Float64(tspan[1]), Float64(tspan[2])), p)
        local sol
        try
            sol = solve(prob, AutoTsit5(Rosenbrock23());
                        callback = landmark_callbacks(log, thresholds;
                                                      cytokinesis = cytokinesis),
                        maxiters = MAXITERS_PER_CELL,
                        save_everystep = false, save_start = false, save_end = true)
        catch
            fates[i] = :Failed
            continue
        end
        if !(string(sol.retcode) in ("Success", "Terminated"))
            fates[i] = :Failed
            continue
        end
        cyc = classify_cycles(log; window = window, cytokinesis = cytokinesis)
        if isempty(cyc)
            fates[i] = :Quiescent
        else
            fates[i] = cyc[1].fate
            d = cyc[1].sg2m_fucci
            d === nothing || (durations[i] = d)
        end
    end

    nfail = count(==(:Failed), fates)
    scored = n - nfail
    counts = Dict(f => count(==(f), fates) for f in labels)
    fracs  = Dict(f => scored > 0 ? counts[f] / scored : NaN for f in labels)
    # Binomial Monte-Carlo error on each fraction.
    mcse   = Dict(f => scored > 0 ? sqrt(max(fracs[f], 0) * (1 - fracs[f]) / scored) : NaN
                  for f in labels)

    d = filter(isfinite, durations)
    dur_cv = length(d) > 1 ? Statistics.std(d) / Statistics.mean(d) : NaN

    return (fractions = fracs, counts = counts, mcse = mcse,
            n = n, scored = scored, failed = nfail,
            duration_cv = dur_cv, duration_mean = isempty(d) ? NaN : Statistics.mean(d),
            sigma = sigma)
end

"""
    calibrate_sigma(; target_cv, sigmas, n, kwargs...) -> NamedTuple

Choose the heterogeneity width by matching the **measured** duration CV.

One number, set by one measurement, and deliberately not by the fate fractions — those
are the hold-out, and fitting the noise to them would make the headline prediction
circular.

Returns the best `sigma`, the CV it achieves, and the full scan so the choice is visible
rather than asserted.
"""
function calibrate_sigma(; target_cv::Real = MEASURED_DURATION_CV,
                           sigmas = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30),
                           n::Int = 400, kwargs...)
    scan = NamedTuple[]
    for s in sigmas
        r = run_ensemble(; n = n, sigma = s, kwargs...)
        push!(scan, (sigma = s, cv = r.duration_cv, scored = r.scored, failed = r.failed))
    end
    ok = filter(x -> isfinite(x.cv), scan)
    isempty(ok) && error("no sigma produced a finite duration CV")
    best = ok[argmin([abs(x.cv - target_cv) for x in ok])]
    return (sigma = best.sigma, cv = best.cv, target = target_cv, scan = scan)
end

"""
    required_n(p; rel=0.1, n_obs=570) -> Int

Ensemble size so the Monte-Carlo error on a fraction `p` is `rel` times the binomial
error the experiment itself carries at `n_obs` cells.

The point of the calculation: there is no value in driving Monte-Carlo error far below
the uncertainty in the measurement being compared against.

**The result does not depend on `p`.** Both errors are binomial in the same `p`, so it
cancels exactly:

    n = p(1-p) / (rel * sqrt(p(1-p)/n_obs))^2 = n_obs / rel^2

That is a more useful fact than it first looks. It means one ensemble size buys the same
error ratio for *every* fate at once — there is no need to size against the rarest one.
At Murganti's 570 cells, `rel = 0.1` gives 57,000 runs and holds for the 1.40 %
polyploidization fraction and the 90.35 % quiescent fraction alike. This is where
TODO.md item 4's "1e4-1e5" comes from.

`p` is retained in the signature because the derivation is about a fraction and dropping
it would obscure why the answer is what it is.
"""
required_n(p::Real; rel::Real = 0.1, n_obs::Int = 570) =
    ceil(Int, p * (1 - p) / (rel * sqrt(p * (1 - p) / n_obs))^2)
