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
                             :ks_E2F7_E2F, :ks_E2F8_E2F, :ks_E2F6,
                             :ks_Ect2_E2F, :ks_AurKB_CDK1, :ks_Anln_E2F, :ks_CSPG,
                             :kf_RhoA, :kf_Midbody,
                             :ks_ATM_ROS, :ks_Ccng1_p53, :kf_CCNB_Ccng1,
                             :fix_p53_massbalance)

"""
New species added by Tier 2, appended after the 63 inherited states.

Order is load-bearing — the inherited RHS destructures positionally, so these must come
last and must never be reordered once trajectories are pinned by test.
"""
const TIER2_SPECIES = ("E2F6", "E2F7", "E2F8",                       # step 2
                       "Ect2", "RhoA", "Centralspindlin",            # step 3
                       "AurKB", "Anillin", "Midbody",
                       "Ccng1")                                      # step 5

"""
Parameter preset that switches the cytokinesis arm on.

Provided as a named constant rather than left to each call site so that "the arm is on"
means one thing everywhere. Requires the E2F split to be on as well — Ect2 is E2F-driven
and E2F8-repressed, so with `ks_E2F8_E2F = 0` the repression arm of it is inert.
"""
const CYTOKINESIS_ON = (ks_Ect2_E2F = 0.60, ks_AurKB_CDK1 = 1.20,
                        ks_Anln_E2F = 0.80, ks_CSPG = 1.00,
                        kf_RhoA = 6.00, kf_Midbody = 4.00)

"""Parameter preset switching on the E2F sub-family split (step 2)."""
const E2F_SPLIT_ON = (ks_E2F7_E2F = 0.20, ks_E2F8_E2F = 0.20)

"""
Parameter preset switching on the maturation axis (step 4) — the **fate** knob.

Only the Ect2 coupling. `maturation_gain = 5.0` puts the abscission arm's closure inside
the observed maturation range rather than outside it; it is a provisional placeholder,
not a fit, and it is the one free scale the axis costs (see `MATURATION_SLOPE_*`).

Deliberately separate from [`E2F6_EXIT_ON`](@ref). The two maturation couplings act on
different observables and calibrate against different data — Ect2 sets the
division/binucleation split against the fate fractions, E2F6 sets cell-cycle exit against
the cycling fractions — so folding them into one preset would tie two independent Phase 3
targets together. Measured with both on at once: everything binucleates and the period
stretches to 100 h, because the knobs fight.
"""
const MATURATION_ON = (maturation_gain = 5.0,)

"""
Parameter preset for E2F6 as the cell-cycle **exit** enforcer (Tier 1's label for it).

Its effect is on the period, and it is strong and graded: at M = 0.55 the cycle runs
39.3 h at `ks_E2F6 = 0`, 48.1 h at 0.01, 57.7 h at 0.02 and 85.6 h at 0.05. That is
maturation-driven cell-cycle exit, the phenomenon the whole project is about, so where it
sits is a Phase 3 calibration against the corrected cycling fractions (P0 17.7 %,
P7 5.5 %) rather than something to choose here. 0.01 is a mild, provisional default.
"""
const E2F6_EXIT_ON = (ks_E2F6 = 0.01,)

"""
Parameter preset switching on the oxidative-stress / DDR arm (step 5).

`ks_p21` is an INHERITED parameter and is raised here from its shipped `1e-4`. That is
not a refit: at `1e-4` p21 is identically zero in every published figure, so the value is
an off-switch rather than an estimate — the source paper's Figure 6 activates the arm the
same way, by hand. The same applies to `kf_ATMp = 0`, which is why ROS drives ATM through
a new term instead.
"""
const DDR_ON = (ks_ATM_ROS = 0.40, ks_Ccng1_p53 = 0.80,
                kf_CCNB_Ccng1 = 200.0, ks_p21 = 0.05,
                fix_p53_massbalance = 1.0)

"""
Every module on, with the fitted parameters at their estimated values.

The only difference from `merge(E2F_SPLIT_ON, CYTOKINESIS_ON, MATURATION_ON, DDR_ON)` is
`ks_Ect2_E2F`, raised 0.60 -> 1.00.

## What that number was fitted to, and what it does not buy

It is fitted to an established fact rather than to a held-out target: neonatal mouse
cardiomyocytes divide around P0-P1 and have largely stopped by P7. With the DDR arm on and
`maturation_gain = 5.0`, the midbody peaks at 0.0579 at P0 and 0.0435 at P7 — a 1.33x
separation straddling the 0.050 abscission threshold. Before this was set, both sat below
it, the model produced binucleation at both ages, and the E2f7/E2f8 knockdown prediction
was vacuous. `ks_Ect2_E2F` in roughly [0.85, 1.15] places the threshold between them; 1.00
is the middle.

**The window is narrow — a ~1.3x range — because the P0/P7 midbody separation is only
1.33x.** That is the same exposure recorded for the abscission threshold in step 3, and it
means the P0/P7 contrast is a threshold-placement result, not a robust one. Two
consequences worth stating rather than discovering later:

  * the knockdown *delta* remains a prediction, but its magnitude is conditioned on where
    the WT sits relative to the threshold, and that position is fitted;
  * profiling the abscission threshold — Tier 1 did this for its own `r080` switch and
    called it "the model's largest single sensitivity" — is now a prerequisite for
    believing any number that comes out of this preset, not an optional extra.
"""
const CALIBRATED = (ks_E2F7_E2F = 0.20, ks_E2F8_E2F = 0.20,
                    ks_Ect2_E2F = 1.00, ks_AurKB_CDK1 = 1.20,
                    ks_Anln_E2F = 0.80, ks_CSPG = 1.00,
                    kf_RhoA = 6.00, kf_Midbody = 4.00,
                    maturation_gain = 5.0,
                    ks_ATM_ROS = 0.40, ks_Ccng1_p53 = 0.80,
                    kf_CCNB_Ccng1 = 200.0, ks_p21 = 0.05,
                    fix_p53_massbalance = 1.0)

"""
    tier2_state(; kwargs...) -> ComponentVector

The 63 inherited states in their original order, plus the E2F sub-family repressors.

Order matters and is load-bearing: the inherited RHS destructures positionally, so the
new components must come last.
"""
function tier2_state(; kwargs...)
    base = NamedTuple(state())
    added = (E2F6 = 0.0, E2F7 = 0.0, E2F8 = 0.0,
             Ect2 = 0.0, RhoA = 0.0, Centralspindlin = 0.0,
             AurKB = 0.0, Anillin = 0.0, Midbody = 0.0,
             Ccng1 = 0.0)
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
        # Repression constants. Raised from an initial 0.10/0.20 after measurement: at
        # those values an E2F6 level of only 0.06 arrests the oscillator outright, which
        # leaves no range to calibrate cell-cycle exit over. Inert while the repressors
        # are absent, so this does not touch the reduction property.
        Ki_E2F78      = 5.00,    # repression constant for the E2F7/E2F8 pool
        Ki_E2F6       = 10.00,   # E2F6 represses more weakly (Tier 1's split)

        # --- step 3: the cytokinesis arm. Enable switches first. ---
        ks_Ect2_E2F   = 0.0,   # Ect2 transcription, E2F-driven
        ks_AurKB_CDK1 = 0.0,   # AurKB (CPC), mitotic
        ks_Anln_E2F   = 0.0,   # anillin transcription
        ks_CSPG       = 0.0,   # centralspindlin bundling
        kf_RhoA       = 0.0,   # Ect2 GEF activity on RhoA
        kf_Midbody    = 0.0,   # furrow -> midbody

        # structural constants, inert while the switches above are zero
        Ki_Ect2_E2F8  = 0.15,  # E2F8 represses Ect2 (Tier 1 r071: E2Fact & !Mat & !E2F8)
        kd_Ect2       = 0.25,
        kd_AurKB      = 0.20,
        # CPC recruitment is switch-like in MPF, not proportional to it -- see the
        # AurKB equation. Ki is well below the unbraked MPF peak (~0.85) so a cell that
        # enters mitosis at all recruits a near-maximal CPC.
        Ki_AurKB_CDK1 = 0.15,
        n_AurKB       = 3.0,
        kd_AurKB_CDH1 = 3.00,  # CPC destroyed at mitotic exit
        kd_Anln       = 0.20,
        kd_Anln_CDH1  = 2.00,
        kd_CSPG       = 0.80,
        Ki_CSPG_CDK1  = 0.05,  # CDK1 phosphorylation blocks centralspindlin bundling
        Ki_RhoA_CDK1  = 0.05,  # and blocks Ect2 GEF activity until anaphase
        kr_RhoA       = 2.50,  # RhoGAP
        RhoA_tot      = 1.00,  # conserved RhoA pool
        kr_Midbody    = 0.60,

        # --- step 4: the maturation axis. See MATURATION_SLOPE_* for the constraint. ---
        M               = 0.0,   # maturation coordinate in [0,1]; 0 => no maturation term
        maturation_gain = 1.0,   # the single free scale for BOTH couplings

        # --- step 5: oxidative stress -> DDR -> Ccng1 -> mitotic-entry brake ---
        ROSenv          = 0.0,   # environmental oxidative stress (input)
        InVitro         = 0.0,   # culture flag (input); Tier 1 raises ROS in culture
        w_InVitro_ROS   = 1.0,   # how much culture adds to ROS
        ks_ATM_ROS      = 0.0,   # enable switch: ROS -> ATM activation
        ks_Ccng1_p53    = 0.0,   # enable switch: p53p -> Ccng1
        kf_CCNB_Ccng1   = 0.0,   # enable switch: Ccng1 -> inhibitory phosphorylation of MPF
        kd_Ccng1        = 0.30,
        # Defaults OFF so the reduction property stays UNCONDITIONAL -- bit-exact for
        # any state, not merely for the dynamically reachable ones where Chk2p = 0.
        # Travels with DDR_ON, which is exactly when the defect can bite.
        fix_p53_massbalance = 0.0,
    )
    return ComponentVector{Float64}(merge(merge(base, added), kwargs))
end

# ---------------------------------------------------------------------------
# The maturation axis, and why it costs ONE parameter rather than two.
#
# Baniol measured two couplings within cycling ventricular cardiomyocytes (n = 89, P0 and
# P7 pooled so developmental stage cannot be doing the work):
#
#     Ect2 ~ M   r = -0.563          E2f6 ~ M   r = +0.396
#
# `../MODEL.md` calls these "the model's two load-bearing couplings", and making them
# measured functions "removes its largest degree of freedom".
#
# But M is defined as mean z(FAO) - mean z(glycolysis), z-scored WITHIN one 285-cell
# dataset, so — as MODEL.md states in its limitations — it "has no absolute cross-system
# scale". A correlation between two z-scored quantities is a regression slope in z-units,
# so what the measurement actually fixes is:
#
#   * the SIGN of each coupling, and
#   * their RATIO, -0.563 / +0.396 = -1.422, which is scale-free and therefore transfers
#     across systems even though neither slope does.
#
# It does NOT fix the absolute strength. So the honest parameterisation is one shared
# `maturation_gain` with the ratio welded in, not two independent slopes. That is one
# fitted parameter for the whole axis instead of two, and it is the difference between
# using the measurement and merely citing it.
# ---------------------------------------------------------------------------

"""Baniol Ect2~M within cycling vCM (n = 89). Frozen: measured, not fitted."""
const MATURATION_SLOPE_ECT2 = -0.563

"""Baniol E2f6~M within cycling vCM (n = 89). Frozen: measured, not fitted."""
const MATURATION_SLOPE_E2F6 = 0.396

"""
    maturation_factors(M, gain) -> (ect2, e2f6)

Multipliers on Ect2 and E2F6 synthesis at maturation `M`.

Both are 1.0 at `M = 0`, which is what makes the axis inert by default. Ect2 is
suppressed through a saturating denominator rather than a linear term, so it cannot go
negative at high `M` — `adult` sits at M = 0.95 and a linear form with any gain above
~1.9 would drive synthesis below zero.
"""
@inline function maturation_factors(M, gain)
    ect2 = 1.0 / (1.0 - gain * MATURATION_SLOPE_ECT2 * M)   # slope is negative -> suppression
    e2f6 = 1.0 + gain * MATURATION_SLOPE_E2F6 * M
    return ect2, e2f6
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
    # E2F6 sits on its own maturation-driven branch, not the activator programme:
    # Baniol find E2f1~E2f6 is only +0.18 n.s., while E2f6~M is +0.396 and E2F6 is the
    # only family member expressed in noncycling cells.
    f_ect2, f_e2f6 = maturation_factors(p.M, p.maturation_gain)
    d.E2F6 = (p.ks_E2F6 * f_e2f6 - p.kd_E2F6 * E2F6) * α

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

    # ------------------------------------------------------------------
    # 6. The cytokinesis arm — Tier 1's central claim, made mechanistic.
    #
    #     Ect2 -> RhoA -> Midbody          (RhoA OBLIGATORY)
    #             Anillin ---^
    #     AurKB -> Centralspindlin --^
    #
    # ## RhoA is obligatory by construction, not by parameter choice
    #
    # `d.Midbody` has exactly ONE production term and RhoA is a factor in it, so no
    # setting of any parameter can make a midbody without RhoA. Tier 1 hit the opposite
    # arrangement four separate times; the worst instance OR'd a second
    # (Centralspindlin & AurKB) route onto the same node, which bypassed the Ect2/RhoA
    # arm entirely and made it *structurally impossible* for Ect2 to be rate-limiting —
    # the model could not express its own central claim. Removing it widened the arm's
    # hiPSC-to-P1 span from 1.6x to 42x. Centralspindlin and AurKB therefore act
    # UPSTREAM, through RhoA, never in parallel with it.
    #
    # ## Timing: why the arm fires after anaphase and not before
    #
    # Ect2's GEF activity and centralspindlin bundling are both blocked by CDK1
    # phosphorylation, so the arm is held off until MPF is destroyed. That is what makes
    # the midbody a late-mitotic structure here rather than something that accumulates
    # through G2, and it is why `Ki_*_CDK1` are small.
    #
    # ## No back-coupling into the core oscillator, and that is a stated limitation
    #
    # The arm reads the cycle but does not drive it: a cell that fails abscission keeps
    # the same molecular oscillation and simply ends up binucleate. That keeps the
    # reduction property trivially exact, but it means ploidy does not feed back on the
    # cycle. Real polyploid cardiomyocytes cycle differently, so this is a Phase 4
    # question, not a permanent claim.
    Ect2, RhoA, CSPG = u.Ect2, u.RhoA, u.Centralspindlin
    AurKB, Anillin, Midbody = u.AurKB, u.Anillin, u.Midbody
    cdk1 = u.CCNB_CDK1

    # Ect2: E2F-driven, repressed by E2F8. Tier 1's r071 is `E2Fact & !Maturation &
    # !E2F8 => Ect2`; the !Maturation arm arrives in step 4. This is where the E2F split
    # earns its keep — the repressor it needs already exists.
    # Tier 1's r071 is `E2Fact & !Maturation & !E2F8 => Ect2`. All three arms are now
    # present: the activator pool, the E2F8 repressor from step 2, and maturation.
    d.Ect2 = (p.ks_Ect2_E2F * E2F * rep * f_ect2 / (1 + E2F8 / p.Ki_Ect2_E2F8)
              - p.kd_Ect2 * Ect2) * α

    # AurKB: chromosomal passenger complex. Mitotic, destroyed by APC/C-Cdh1 at exit.
    #
    # Recruitment is SWITCH-LIKE in MPF, not proportional to it. Driving it linearly was
    # a modelling error with a measurable consequence: the Ccng1 brake lowers MPF, so a
    # linear AurKB drive propagated that reduction straight down
    # AurKB -> Centralspindlin -> RhoA -> Midbody and the DDR arm silently suppressed
    # CYTOKINESIS as a side effect of braking mitotic ENTRY. Division then became
    # unreachable in every mouse context and the E2f7/E2f8 knockdown prediction went
    # vacuous (dDivision = 0.0000 at both P0 and P7).
    #
    # The biology says otherwise: the CPC is recruited to centromeres at mitotic entry
    # and stays active through mitosis largely independently of the precise MPF level. A
    # Hill term keeps the property that matters -- no mitosis, no CPC, hence no midbody
    # (factor 0.04 at MPF 0.05) -- while decoupling a cell that DOES enter mitosis with
    # reduced MPF from losing its midbody (0.73 at MPF 0.21, against 0.21 linear).
    cpc = cdk1^p.n_AurKB / (p.Ki_AurKB_CDK1^p.n_AurKB + cdk1^p.n_AurKB)
    d.AurKB = (p.ks_AurKB_CDK1 * cpc
               - p.kd_AurKB * AurKB
               - p.kd_AurKB_CDH1 * AurKB * u.APCC_CDH1) * α

    # Centralspindlin: bundled by AurKB, blocked by CDK1 until anaphase.
    d.Centralspindlin = (p.ks_CSPG * AurKB / (1 + cdk1 / p.Ki_CSPG_CDK1)
                         - p.kd_CSPG * CSPG) * α

    # Anillin: furrow scaffold, E2F-driven, cleared at mitotic exit.
    d.Anillin = (p.ks_Anln_E2F * E2F * rep
                 - p.kd_Anln * Anillin
                 - p.kd_Anln_CDH1 * Anillin * u.APCC_CDH1) * α

    # RhoA-GTP: activated by Ect2 at the central spindle once CDK1 has fallen.
    # Written against a conserved pool so activation saturates.
    d.RhoA = (p.kf_RhoA * Ect2 * CSPG / (1 + cdk1 / p.Ki_RhoA_CDK1)
              * (p.RhoA_tot - RhoA)
              - p.kr_RhoA * RhoA) * α

    # Midbody. One production term; RhoA and anillin both required.
    d.Midbody = (p.kf_Midbody * RhoA * Anillin - p.kr_Midbody * Midbody) * α

    # ------------------------------------------------------------------
    # 7. Oxidative stress -> DDR -> Ccng1 -> the mitotic-entry brake.
    #
    # Tier 1's other half of the 2x2: "maturation closes the abscission arm, and culture
    # closes the mitotic-entry brake (ROS -> DDR -> Ccng1/Pkmyt1 -> MitoticEntry,
    # collapsing mitotic entry from 0.531 in vivo to 0.251 in vitro). Which of the two
    # shuts first is the whole 2x2." Step 4 built the maturation half; this is the rest.
    #
    # ## Turning the arm on is not a refit of a published parameter
    #
    # The inherited model SHIPS this arm dormant: `kf_ATMp = 0`, with the damage input
    # commented out in `d.ATM` as `#-kf_ATMp*(KDDS+sig)`, and `ks_p21 = 1e-4`, which
    # leaves p21 identically zero in every published figure. Those are off-switches, not
    # fitted values -- Figure 6 of the source paper activates the arm by hand the same
    # way. Raising them is enabling shipped-but-unused machinery, not re-estimating
    # something the manuscript reports.
    #
    # ROS enters as an input rather than a state: it is an environmental property of the
    # culture, constant on the timescale of a cycle, exactly like M.
    ros = p.ROSenv + p.w_InVitro_ROS * p.InVitro
    if p.ks_ATM_ROS != 0.0 && ros != 0.0
        atm_act = p.ks_ATM_ROS * ros * u.ATM * α
        d.ATM  -= atm_act          # mass-conserving: ATM -> ATMp
        d.ATMp += atm_act
    end

    # p53 <-> p53p mass balance. `d.p53p` in the inherited model reads
    #     (kf_p53p*Chk2p) - kr_p53p*PPase*p53p
    # while `d.p53`'s matching loss term is `-(kf_p53p*Chk2p)*p53`. The `*p53` factor is
    # missing on the gain side, so phosphorylation adds p53p at a rate independent of how
    # much p53 there is, and total p53 is not conserved. Harmless in the published
    # figures because Chk2p is identically zero there -- which is also why applying the
    # fix unconditionally does not disturb the reduction property, asserted by test.
    if p.fix_p53_massbalance != 0.0
        d.p53p += p.fix_p53_massbalance * p.kf_p53p * u.Chk2p * (u.p53 - 1.0) * α
    end

    # Ccng1: p53-induced, and the brake on mitotic entry. Implemented as an extra
    # inhibitory phosphorylation of MPF -- written as a transfer between CCNB_CDK1 and
    # CCNB_CDK1p so it cannot create or destroy MPF, only inactivate it.
    d.Ccng1 = (p.ks_Ccng1_p53 * u.p53p - p.kd_Ccng1 * u.Ccng1) * α
    if p.kf_CCNB_Ccng1 != 0.0
        brake = p.kf_CCNB_Ccng1 * u.Ccng1 * cdk1 * α
        d.CCNB_CDK1  -= brake
        d.CCNB_CDK1p += brake
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
                       cytokinesis::Bool = false,
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
                    callback = landmark_callbacks(log, thresholds;
                                                  cytokinesis = cytokinesis), kwargs...)
        return sol, log
    end
    return solve(prob, AutoTsit5(Rosenbrock23()); kwargs...)
end
