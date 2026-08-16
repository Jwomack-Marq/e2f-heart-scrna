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
end
EventLog() = EventLog(Float64[], Float64[], Float64[], Float64[], Float64[], Float64[])

Base.isempty(l::EventLog) = all(isempty, (l.restriction, l.s_entry, l.neb,
                                          l.anaphase, l.mitotic_exit, l.geminin_on))

function Base.show(io::IO, l::EventLog)
    print(io, "EventLog(restriction=", length(l.restriction),
              ", s_entry=", length(l.s_entry),
              ", neb=", length(l.neb),
              ", anaphase=", length(l.anaphase),
              ", mitotic_exit=", length(l.mitotic_exit),
              ", geminin_on=", length(l.geminin_on), ")")
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
"""
function landmark_callbacks(log::EventLog, thr::EventThresholds = EventThresholds())
    up(idx, level, sink) = ContinuousCallback(
        (u, t, integ) -> u[idx] - level,
        integ -> push!(sink, integ.t),   # up-crossing
        nothing)
    down(idx, level, sink) = ContinuousCallback(
        (u, t, integ) -> u[idx] - level,
        nothing,
        integ -> push!(sink, integ.t))   # down-crossing

    # Total geminin is a sum of two states, so it needs its own condition rather than a
    # single index. This is the green FUCCI channel appearing — the landmark Murganti
    # Fig 2E actually times from.
    ig, igc = species_index("Geminin"), species_index("Geminin_CDT1")
    geminin_cb = ContinuousCallback(
        (u, t, integ) -> u[ig] + u[igc] - thr.geminin_on,
        integ -> push!(log.geminin_on, integ.t),
        nothing)

    return CallbackSet(
        up(species_index("ppRB"),      thr.restriction,  log.restriction),
        up(species_index("CCNE_CDK2"), thr.s_entry,      log.s_entry),
        up(species_index("LMNAp"),     thr.neb,          log.neb),
        down(species_index("PTTG1"),   thr.anaphase,     log.anaphase),
        down(species_index("CCNB_CDK1"), thr.mitotic_exit, log.mitotic_exit),
        geminin_cb,
    )
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
                             kwargs...)
    log = EventLog()
    u0 = state()
    p = copy(params())
    p.α = alpha
    p.con_VOL = con_VOL
    p.con_ABE = con_ABE
    prob = ODEProblem(modelDiffEq!, u0, (Float64(tspan[1]), Float64(tspan[2])), p)
    sol = solve(prob, AutoTsit5(Rosenbrock23());
                callback = landmark_callbacks(log, thresholds), kwargs...)
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
                    keep(log.anaphase), keep(log.mitotic_exit), keep(log.geminin_on))
end
