# The Tier-2 extended model: E2F sub-family split, and FUCCI licensing coupled to E2F.
#
# ## How this preserves the published model exactly
#
# The inherited `modelDiffEq!` is **called verbatim, never copied**. `diff_eqns.jl`
# destructures `u` and `p` positionally, so an extended `ComponentVector` whose extra
# components are appended at the END is destructured to exactly the same 63 states and
# 218 parameters — verified bit-identical, and asserted by test.
#
# Tier 2 is then a set of ADDITIVE CORRECTIONS on top of that call, each written so it is
# identically zero at default parameters:
#
#     d.CCNE += ks_CCNE_E2F * E2F * (rep - 1)        # zero when rep == 1
#     d.CDT1 += ks_CDT1_E2F * w * (E2F/E2F_ref - 1)  # zero when w == 0
#
# So `tier2DiffEq!` reduces to `modelDiffEq!` **structurally**, not merely to within a
# tolerance. Copying the 442-line RHS and editing it would have recreated exactly the
# hazard flagged in Phase 0 — `diff_eqn_plk1p_test.jl`, a shadow copy that
# `include_all_scripts()` could let silently win.
#
# Every correction is multiplied by `α`, as the inherited RHS multiplies every term, so
# enabling a module cannot silently rescale the clock.

"""
    TIER2_ENABLE_PARAMS

The five parameters that gate every Tier-2 addition. All zero (the default) means the
extended model *is* the published model. Nothing else needs to be switched: the
remaining new parameters are rate constants that multiply species which stay at zero
while these are zero, so they can carry biologically motivated defaults without
affecting the reduction.
"""
const TIER2_ENABLE_PARAMS = (:w_CDT1_E2F, :w_Geminin_E2F, :kd_CDT1_CDK2,
                             :ks_E2F7_E2F, :ks_E2F8_E2F, :ks_E2F6)

"""New species added by Tier 2, appended after the 63 inherited states."""
const TIER2_SPECIES = ("E2F6", "E2F7", "E2F8")

"""
    tier2_state(; kwargs...) -> ComponentVector

The 63 inherited states in their original order, plus the E2F sub-family repressors.

Order matters and is load-bearing: the inherited RHS destructures positionally, so the
new components must come last.
"""
function tier2_state(; kwargs...)
    base = NamedTuple(state())
    added = (E2F6 = 0.0, E2F7 = 0.0, E2F8 = 0.0)
    return ComponentVector{Float64}(merge(merge(base, added), kwargs))
end

"""All state names, inherited then Tier 2. The first 63 match [`state_names`](@ref)."""
tier2_state_names() = vcat(state_names(), collect(TIER2_SPECIES))

"""
    tier2_params(; kwargs...) -> ComponentVector

The 218 published parameters, unchanged and in their original order, plus the Tier-2
block appended.

## Provenance of the new values

| parameter | default | where it comes from |
|---|---|---|
| `w_CDT1_E2F`, `w_Geminin_E2F` | 0 | **enable switches** — fraction of synthesis that is E2F-driven |
| `ks_E2F7_E2F`, `ks_E2F8_E2F` | 0 | **enable switches** — E2F1-driven induction |
| `ks_E2F6` | 0 | **enable switch** — maturation-driven in step 4, constant until then |
| `E2F_ref` | 0.1218 | measured: cycle-mean free E2F at the published α, so total licensing synthesis is conserved at `w = 1` rather than silently rescaled |
| `kd_E2F7`, `kd_E2F8`, `kd_E2F6` | 0.30 / 0.90 / 0.10 | E2F8 turns over fastest: Baniol's `E2f8~Ccna2` = −0.56 makes it G1/S-restricted. E2F6 is the stable one — the only family member expressed in noncycling cells (0.38 in 100 % of noncycling P7 vCM) |
| `kd_E2F8_CCNA` | 2.0 | the explicit G1/S restriction: CycA-CDK2 clears E2F8 as cells pass into S/G2 |
| `Ki_E2F78`, `Ki_E2F6` | 0.10 / 0.20 | repression constants; E2F6 is the weaker repressor per Tier 1's split |

Nothing here is fitted yet. These are structural defaults chosen so the module behaves
before calibration; the fitted budget is Phase 3's problem and is declared there.
"""
function tier2_params(; kwargs...)
    base = NamedTuple(params())
    added = (
        # --- enable switches: all zero => exact reduction to the published model ---
        w_CDT1_E2F    = 0.0,   # fraction of Cdt1 synthesis driven by E2F
        w_Geminin_E2F = 0.0,   # fraction of geminin synthesis driven by E2F
        kd_CDT1_CDK2  = 0.0,   # CDK2-driven Cdt1 clearance at S onset (see below)
        ks_E2F7_E2F   = 0.0,   # E2F1 -> E2F7 induction (Baniol r = +0.49)
        ks_E2F8_E2F   = 0.0,   # E2F1 -> E2F8 induction (Baniol r = +0.46)
        ks_E2F6       = 0.0,   # E2F6 synthesis; maturation-driven from step 4

        # --- structural constants, inert while the switches are zero ---
        E2F_ref       = 0.1218,  # cycle-mean free E2F, measured at the published alpha
        kd_E2F7       = 0.30,
        kd_E2F8       = 0.90,
        kd_E2F6       = 0.10,
        kd_E2F8_CCNA  = 2.00,    # CycA-CDK2 clears E2F8 -> G1/S restriction
        Ki_E2F78      = 0.10,    # repression constant for the E2F7/E2F8 pool
        Ki_E2F6       = 0.20,    # E2F6 represses more weakly
    )
    return ComponentVector{Float64}(merge(merge(base, added), kwargs))
end

"""
    e2f_repression(E2F6, E2F7, E2F8, p) -> Float64

Fractional repression of E2F-driven transcription by the repressor sub-family, in
`(0, 1]`. Exactly 1 when all three repressors are absent.

E2F7 and E2F8 act as one pool because they are obligate homo/heterodimers on the same
sites and Tier 1 found an OR'd pool of repressors saturates; E2F6 is separate because
Baniol place it on its own maturation-driven branch rather than the activator programme
(it is the only member expressed in noncycling cells, and `E2f1~E2f6` is only +0.18 n.s.).

A single saturating denominator, not a product of independent factors — competitive
occupancy of the same promoters, so two repressors do not multiply into an
implausibly deep block.
"""
@inline function e2f_repression(E2F6, E2F7, E2F8, Ki_E2F78, Ki_E2F6)
    return 1.0 / (1.0 + (E2F7 + E2F8) / Ki_E2F78 + E2F6 / Ki_E2F6)
end

"""
    tier2DiffEq!(d, u, p, t)

Extended right-hand side: the published model plus Tier-2 corrections.

Reduces to `modelDiffEq!` identically when the five [`TIER2_ENABLE_PARAMS`](@ref) are
zero — structurally, because every correction carries a factor that is then zero.
"""
function tier2DiffEq!(d, u, p, t)
    # 1. The published model, verbatim, on the leading 63 states / 218 parameters.
    modelDiffEq!(d, u, p, t)

    α = p.α
    E2F = u.E2F
    E2F6, E2F7, E2F8 = u.E2F6, u.E2F7, u.E2F8

    # 2. The E2F sub-family repressors. Inert while their synthesis rates are zero.
    #    E2F8 additionally carries CycA-driven clearance, which is what makes it
    #    G1/S-restricted rather than merely short-lived.
    d.E2F7 = (p.ks_E2F7_E2F * E2F - p.kd_E2F7 * E2F7) * α
    d.E2F8 = (p.ks_E2F8_E2F * E2F - p.kd_E2F8 * E2F8
              - p.kd_E2F8_CCNA * E2F8 * u.CCNA_CDK2) * α
    d.E2F6 = (p.ks_E2F6 - p.kd_E2F6 * E2F6) * α

    # 3. Repression of E2F-driven transcription. `rep - 1` is zero when unrepressed, so
    #    each of these vanishes at default parameters.
    rep = e2f_repression(E2F6, E2F7, E2F8, p.Ki_E2F78, p.Ki_E2F6)
    Δrep = rep - 1.0
    if Δrep != 0.0
        # The six consumers of free E2F, matching diff_eqns.jl lines 131/165/175/245/322/335.
        # `ks_CCNE_pRBE2F * pRB_E2F` is deliberately NOT repressed: that is E2F held on
        # phospho-Rb, not the free pool the repressors compete with.
        d.CCNE    += p.ks_CCNE_E2F   * E2F * Δrep * α
        d.CDC25A  += p.ks_CDC25A_E2F * E2F * Δrep * α
        d.CCNA    += p.ks_CCNA_E2F   * E2F * Δrep * α
        d.PLK1    += p.ks_PLK1 * E2F * (1 / (1 + u.p53p / p.ki_PLK1)) * Δrep * α
        d.PTTG1   += p.ks_PTTG1      * E2F * Δrep * α
        d.EMI1    += p.ks_EMI1_E2F   * E2F * Δrep * α
        # E2F auto-repression: the canonical E2F7/8 negative feedback loop. Applies to
        # the basal synthesis term, which is what `d.E2F` actually carries.
        d.E2F     += p.ks_E2F * Δrep * α
    end

    # 4. Licensing coupled to E2F.
    #    Phase 1 measured the defect this addresses: total geminin clears the published
    #    0.05 FUCCI cutoff for only 12.6 % of the cycle (mAG should mark ~40 %), and Cdt1
    #    and geminin never overlap, so the G1/S double-positive has frequency zero
    #    against 1.6 % (Murganti Fig 1C) and 19.1 % (Baniol Fig 1D at P0).
    #
    #    Cause: ks_CDT1_E2F and ks_Geminin_E2F are CONSTANTS in diff_eqns.jl:425,432
    #    despite their names and comments, so licensing is decoupled from the cycle.
    #    Here synthesis becomes (1-w) basal + w E2F-driven, normalised by E2F_ref so
    #    total synthesis is conserved at w = 1 instead of being silently rescaled.
    #    Repression applies to the E2F-driven part only.
    if p.w_CDT1_E2F != 0.0
        d.CDT1    += p.ks_CDT1_E2F * p.w_CDT1_E2F * (E2F / p.E2F_ref * rep - 1.0) * α
    end
    if p.w_Geminin_E2F != 0.0
        d.Geminin += p.ks_Geminin_E2F * p.w_Geminin_E2F * (E2F / p.E2F_ref * rep - 1.0) * α
    end

    # 5. CDK2-driven Cdt1 clearance at S-phase onset.
    #    Coupling licensing to E2F is necessary but not sufficient, and the measurement
    #    says why: with only steps 1-2 the geminin window sits ENTIRELY INSIDE the Cdt1
    #    window (mAG-only frequency 0.000), so there is a G1/S double-positive but no
    #    S/G2/M state at all. Cdt1's only degradation route in the inherited model is
    #    SCF, and SCF peaks at mitosis rather than at the G1/S transition, so nothing
    #    clears Cdt1 when S begins.
    #
    #    The missing edge is real and well documented: CDK2 phosphorylation of Cdt1
    #    (T29) licenses SCF-Skp2 recognition, and CRL4-Cdt2 degrades PCNA-loaded Cdt1
    #    through S phase. Both are S-phase-restricted, which is exactly the timing the
    #    reporter needs. Written against CycE-CDK2 + CycA-CDK2, the two active CDK2
    #    pools, in the mass-action style of the surrounding equations.
    if p.kd_CDT1_CDK2 != 0.0
        d.CDT1 -= p.kd_CDT1_CDK2 * u.CDT1 * (u.CCNE_CDK2 + u.CCNA_CDK2) * α
    end

    return nothing
end

"""
    solve_tier2(; alpha, tspan, enable, thresholds, record_events, kwargs...)

Solve the extended model. `enable` is a NamedTuple of Tier-2 parameter overrides, e.g.
`(w_Geminin_E2F = 1.0, ks_E2F7_E2F = 0.5)`. With no overrides this is the published
model, and the returned trajectory matches [`solve_baseline`](@ref) exactly on the 63
inherited states.

Returns `sol`, or `(sol, log)` when `record_events = true`.
"""
function solve_tier2(; alpha::Real = PUBLISHED_ALPHA,
                       tspan::Tuple{<:Real,<:Real} = (0.0, 2500.0),
                       enable::NamedTuple = NamedTuple(),
                       con_VOL::Real = 0.0, con_ABE::Real = 0.0,
                       record_events::Bool = false,
                       thresholds::EventThresholds = EventThresholds(),
                       kwargs...)
    u0 = tier2_state()
    p = tier2_params(; enable...)
    p.α = alpha
    p.con_VOL = con_VOL
    p.con_ABE = con_ABE
    prob = ODEProblem(tier2DiffEq!, u0, (Float64(tspan[1]), Float64(tspan[2])), p)
    if record_events
        log = EventLog()
        sol = solve(prob, AutoTsit5(Rosenbrock23());
                    callback = landmark_callbacks(log, thresholds), kwargs...)
        return sol, log
    end
    return solve(prob, AutoTsit5(Rosenbrock23()); kwargs...)
end
