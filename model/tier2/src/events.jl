# Discrete cell-cycle landmarks, detected by root-finding rather than by scanning
# saved points.
#
# `../TODO.md` item 3 asks for "discrete envelope-breakdown and abscission events". This
# is the envelope-breakdown half; abscission needs the Ect2/RhoA/midbody arm and lands
# in Phase 2.
#
# Why callbacks and not post-hoc peak-finding: phase *durations* are the Phase 1 gate
# (the five `duration` rows in ../cmcycle/data/cmcycle_targets.csv), so event times are
# the measurement, not a diagnostic. `ContinuousCallback` root-finds the crossing to
# solver tolerance; scanning the saved grid would quantise every duration to the step
# size, which the adaptive solver varies by orders of magnitude across a cycle.

"""
    EventThresholds

Absolute activity thresholds defining each landmark.

**Deliberately absolute, not relative to each trajectory's own amplitude.** A fraction-of-
own-maximum rule would make every trajectory cross every threshold by construction,
including a near-quiescent cardiomyocyte whose CycE never leaves the noise floor — and
`Quiescent` is one of the four fates that has to be *detectable*, not assumed away.
Holding them fixed across contexts is what lets a cell genuinely fail to enter S.

Defaults are half of each species' peak in the published proliferating baseline
(α = 1.447), except `mitotic_exit`, which is set low because mitotic exit is MPF
*destruction* rather than a midpoint crossing. Derived by `scripts/derive_thresholds.jl`
and pinned by test.
"""
Base.@kwdef struct EventThresholds
    s_entry::Float64      = 0.299   # CCNE_CDK2, 0.50 x baseline peak 0.5970
    neb::Float64          = 0.180   # LMNAp,     0.50 x baseline peak 0.3600
    anaphase::Float64     = 0.165   # PTTG1,     0.50 x baseline peak 0.3298
    mitotic_exit::Float64 = 0.099   # CCNB_CDK1, 0.15 x baseline peak 0.6622
    restriction::Float64  = 0.073   # ppRB,      0.50 x baseline peak 0.1466
    geminin_on::Float64   = FUCCI_THRESHOLD  # total geminin, the published FUCCI cutoff
    abscission::Float64   = 0.050  # Midbody; only used when the cytokinesis arm is on

    # Minimum gap between two recordings of the SAME landmark, in hours.
    #
    # Not cosmetic. A ContinuousCallback on a level that a slowly-varying signal happens
    # to sit near will chatter: with E2F6 repression on, total geminin hovers at the
    # FUCCI cutoff and logged 995,356 crossings in one run, driving the solver into
    # maxiters at t = 1518. Every landmark here is once-per-cycle and the shortest cycle
    # observed is ~28 h, so 1 h is far below any real inter-event gap and far above the
    # chatter.
    refractory::Float64   = 1.0

    # Schmitt-trigger deadband, as a fraction of each level.
    #
    # The refractory guard above cleans the LOG but not the solver: a ContinuousCallback
    # still root-finds and restarts the step at every crossing it detects, so a signal
    # lingering near its level costs real work even when the extra events are discarded.
    # Total geminin spends 12-19 % of each cycle within +/-10 % of the FUCCI cutoff, and
    # that was enough to drive 2e6 steps and a MaxIters abort at t = 1518.
    #
    # With hysteresis the detector must see the signal move away by this fraction before
    # it will fire again, which removes the chatter at source rather than after the fact.
    hysteresis::Float64   = 0.25
end

"""
    EventLog

Times at which each landmark was crossed, in simulation hours.

`s_entry` and `restriction` are up-crossings; `anaphase` and `mitotic_exit` are
down-crossings (securin destruction and MPF destruction). `neb` is an up-crossing of
lamin phosphorylation.
"""
struct EventLog
    restriction::Vector{Float64}
    s_entry::Vector{Float64}
    neb::Vector{Float64}
    anaphase::Vector{Float64}
    mitotic_exit::Vector{Float64}
    geminin_on::Vector{Float64}
    abscission::Vector{Float64}
end
EventLog() = EventLog(Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
                      Float64[])

Base.isempty(l::EventLog) = all(isempty, (l.restriction, l.s_entry, l.neb, l.anaphase,
                                          l.mitotic_exit, l.geminin_on, l.abscission))

function Base.show(io::IO, l::EventLog)
    print(io, "EventLog(restriction=", length(l.restriction),
              ", s_entry=", length(l.s_entry),
              ", neb=", length(l.neb),
              ", anaphase=", length(l.anaphase),
              ", mitotic_exit=", length(l.mitotic_exit),
              ", geminin_on=", length(l.geminin_on),
              ", abscission=", length(l.abscission), ")")
end

"""
    landmark_callbacks(log, thr) -> CallbackSet

Build the callback set that fills `log`.

Each landmark is its own `ContinuousCallback` with the unused direction passed as
`nothing`, so a species that oscillates around a threshold cannot log the same landmark
twice per cycle from the two crossing directions.

The affects record only — they never touch `u` or `p`. Bookkeeping is derived from the
log afterwards (see `fates.jl`), which keeps the ODE trajectory identical to the
callback-free solve. That is asserted by test.

The `cytokinesis` flag adds the abscission landmark. It is off by default and must stay
that way: `Midbody` is a Tier-2 species at an index past the end of the inherited
63-component state, so arming that callback against a `solve_baseline` problem would
index out of bounds.
"""
function landmark_callbacks(log::EventLog, thr::EventThresholds = EventThresholds();
                            cytokinesis::Bool = false)
    # Refractory-guarded recording: drop a repeat of the same landmark that lands within
    # `thr.refractory` hours of the previous one. See EventThresholds.refractory.
    rec!(sink, t) = (isempty(sink) || t - last(sink) > thr.refractory) && push!(sink, t)

    # Schmitt-triggered detectors. `armed` flips the level the condition watches, so once
    # a landmark fires the signal must retreat past the deadband before it can fire
    # again. Safe here because the affects only record -- they never touch u or p -- so
    # switching the condition mid-solve introduces no state discontinuity.
    function up(idx, level, sink; f = (u, i) -> u[i])
        lo = level * (1 - thr.hysteresis)
        armed = Ref(true)
        ContinuousCallback(
            (u, t, integ) -> f(u, idx) - (armed[] ? level : lo),
            integ -> (armed[] && (rec!(sink, integ.t); armed[] = false); nothing),
            integ -> (armed[] = true; nothing))
    end
    function down(idx, level, sink)
        hi = level * (1 + thr.hysteresis)
        armed = Ref(true)
        ContinuousCallback(
            (u, t, integ) -> u[idx] - (armed[] ? level : hi),
            integ -> (armed[] = true; nothing),
            integ -> (armed[] && (rec!(sink, integ.t); armed[] = false); nothing))
    end

    # Total geminin is a sum of two states, so it needs its own condition rather than a
    # single index. This is the green FUCCI channel appearing — the landmark Murganti
    # Fig 2E actually times from.
    ig, igc = species_index("Geminin"), species_index("Geminin_CDT1")
    geminin_cb = up(ig, thr.geminin_on, log.geminin_on;
                    f = (u, i) -> u[i] + u[igc])

    core = (
        up(species_index("ppRB"),      thr.restriction,  log.restriction),
        up(species_index("CCNE_CDK2"), thr.s_entry,      log.s_entry),
        up(species_index("LMNAp"),     thr.neb,          log.neb),
        down(species_index("PTTG1"),   thr.anaphase,     log.anaphase),
        down(species_index("CCNB_CDK1"), thr.mitotic_exit, log.mitotic_exit),
        geminin_cb,
    )
    cytokinesis || return CallbackSet(core...)
    return CallbackSet(core...,
        up(species_index("Midbody"), thr.abscission, log.abscission))
end

"""
    solve_with_events(; alpha, tspan, thresholds, kwargs...) -> (sol, log)

Solve the baseline model while recording cell-cycle landmarks.

Extra keyword arguments are forwarded to `solve`, so a caller can raise `saveat` density
for plotting without affecting the event times — those come from root-finding, not from
the saved grid.
"""
function solve_with_events(; alpha::Real = PUBLISHED_ALPHA,
                             tspan::Tuple{<:Real,<:Real} = (0.0, 2500.0),
                             thresholds::EventThresholds = EventThresholds(),
                             con_VOL::Real = 0.0, con_ABE::Real = 0.0,
                             cytokinesis::Bool = false,
                             kwargs...)
    log = EventLog()
    u0 = state()
    p = copy(params())
    p.α = alpha
    p.con_VOL = con_VOL
    p.con_ABE = con_ABE
    prob = ODEProblem(modelDiffEq!, u0, (Float64(tspan[1]), Float64(tspan[2])), p)
    sol = solve(prob, AutoTsit5(Rosenbrock23());
                callback = landmark_callbacks(log, thresholds;
                                              cytokinesis = cytokinesis), kwargs...)
    return sol, log
end

"""
    trim(log, window) -> EventLog

Restrict every landmark series to `window`, for analysing the settled limit cycle
rather than the start-up transient.
"""
function trim(log::EventLog, window::Tuple{<:Real,<:Real})
    keep(v) = filter(t -> window[1] <= t <= window[2], v)
    return EventLog(keep(log.restriction), keep(log.s_entry), keep(log.neb),
                    keep(log.anaphase), keep(log.mitotic_exit), keep(log.geminin_on),
                    keep(log.abscission))
end
