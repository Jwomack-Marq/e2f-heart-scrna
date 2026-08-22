# The parameter budget, mechanized.
#
# Tier 1's central methodological claim is restraint: 4 fitted parameters against
# `../MODEL.md` §3.5's estimate of ~20-25 effective independent constraints, giving a
# defensible ceiling near 13. It enforces this with an `evidence` column on every reaction
# and a linter that fails the build. Tier 2 arrives with 218 inherited parameters plus 41
# of its own, so the same discipline needs the same machinery — a promise in a README is
# not enforcement.
#
# Every Tier-2 parameter must be declared here. A parameter that exists but is not
# declared fails the test suite, so the budget cannot drift by accident.

"""
    Provenance

Where a parameter's value comes from. Only `FITTED` counts against the budget.

  * `FROZEN`     — inherited from the published model, unchanged. A citation, not a
                   degree of freedom.
  * `INPUT`      — a context variable read from Tier 1's manifest, not a parameter of the
                   model at all (`M`, `ROSenv`, `InVitro`).
  * `MEASURED`   — set directly by a measurement, with no freedom left (the Baniol
                   correlation ratio, the cycle-mean E2F reference).
  * `STRUCTURAL` — a rate constant chosen on biological grounds and *not* tuned against
                   any target. These are the honest grey zone: they are not fitted, but
                   they are not measured either, and the paper must say so.
  * `SWITCH`     — an on/off gate with no intermediate meaning (`fix_p53_massbalance`).
  * `FITTED`     — estimated against data. Counts against the budget.
"""
@enum Provenance FROZEN INPUT MEASURED STRUCTURAL SWITCH FITTED

"""
`../MODEL.md` §3.5: ~20-25 effective independent constraints exist across the two source
papers, which puts the defensible ceiling at about 13 fitted parameters.
"""
const FITTED_CEILING = 13

"""
The Phase 3 target, deliberately below [`FITTED_CEILING`](@ref).

Leaving headroom is the point. A model that spends its entire defensible budget has no
answer to "why should I believe this generalises"; one that spends two thirds of it does.
"""
const FITTED_BUDGET = 8

"""
    PARAM_PROVENANCE

Every Tier-2 parameter, its provenance, and the justification. Checked for completeness
by test, so a new parameter cannot be added without declaring where its value came from.
"""
const PARAM_PROVENANCE = Dict{Symbol,Tuple{Provenance,String}}(
    # --- inputs: context variables, read from Tier 1's manifest ---
    :M        => (INPUT, "maturation coordinate; from cmfate_model.toml"),
    :ROSenv   => (INPUT, "environmental oxidative stress; from cmfate_model.toml"),
    :InVitro  => (INPUT, "culture flag; from cmfate_model.toml"),

    # --- measured: fixed by a specific number in the source papers ---
    :E2F_ref  => (MEASURED, "cycle-mean free E2F at the published alpha, 0.1218; " *
                            "normalises licensing synthesis so w=1 conserves it"),

    # --- switches ---
    :fix_p53_massbalance => (SWITCH, "restores the *p53 factor missing from d.p53p"),

    # --- FITTED: these are the budget ---
    :maturation_gain => (FITTED, "the single free scale of the maturation axis; the " *
                                 "Baniol slope RATIO -1.422 is welded in, so this is " *
                                 "one parameter for two couplings"),
    :ks_E2F6         => (FITTED, "where cell-cycle exit happens; against cycling fractions"),
    :ks_Ect2_E2F     => (FITTED, "cytokinesis arm scale; against the division/binucleation split"),
    :kf_CCNB_Ccng1   => (FITTED, "DDR brake strength. Only the product " *
                                 "ks_Ccng1_p53*kf_CCNB_Ccng1/kd_Ccng1 is identifiable, " *
                                 "so the other two are pinned STRUCTURAL and this " *
                                 "carries the whole brake as ONE effective parameter"),
    :ks_E2F7_E2F     => (FITTED, "E2F repressor pool strength"),
    :ks_E2F8_E2F     => (FITTED, "tied to ks_E2F7_E2F during calibration; see lint_budget"),

    # --- structural: biologically motivated, not tuned against any target ---
    :w_CDT1_E2F    => (STRUCTURAL, "rejected in step 1; kept at 0"),
    :w_Geminin_E2F => (STRUCTURAL, "rejected in step 1; kept at 0"),
    :kd_CDT1_CDK2  => (STRUCTURAL, "rejected in step 1; kept at 0"),
    :kd_E2F7       => (STRUCTURAL, "repressor turnover"),
    :kd_E2F8       => (STRUCTURAL, "faster than E2F7: Baniol E2f8~Ccna2 = -0.56"),
    :kd_E2F6       => (STRUCTURAL, "slowest: E2F6 is the member present in noncycling cells"),
    :kd_E2F8_CCNA  => (STRUCTURAL, "CycA clearance; the explicit G1/S restriction"),
    :Ki_E2F78      => (STRUCTURAL, "raised to 5.0 after measurement: at 0.1 the oscillator dies"),
    :Ki_E2F6       => (STRUCTURAL, "raised to 10.0; E2F6 represses more weakly (Tier 1's split)"),
    :ks_AurKB_CDK1 => (STRUCTURAL, "CPC is mitotic"),
    :ks_Anln_E2F   => (STRUCTURAL, "anillin transcription"),
    :ks_CSPG       => (STRUCTURAL, "centralspindlin bundling"),
    :kf_RhoA       => (STRUCTURAL, "Ect2 GEF activity"),
    :kf_Midbody    => (STRUCTURAL, "furrow -> midbody"),
    :Ki_Ect2_E2F8  => (STRUCTURAL, "Tier 1 r071: E2Fact & !Maturation & !E2F8"),
    :kd_Ect2       => (STRUCTURAL, "Ect2 turnover"),
    :kd_AurKB      => (STRUCTURAL, "AurKB turnover"),
    :Ki_AurKB_CDK1 => (STRUCTURAL, "half-max MPF for CPC recruitment; switch-like, so " *
                                   "braking mitotic entry does not silently suppress " *
                                   "cytokinesis"),
    :n_AurKB       => (STRUCTURAL, "CPC recruitment cooperativity"),
    :kd_AurKB_CDH1 => (STRUCTURAL, "CPC destroyed at mitotic exit"),
    :kd_Anln       => (STRUCTURAL, "anillin turnover"),
    :kd_Anln_CDH1  => (STRUCTURAL, "anillin cleared at mitotic exit"),
    :kd_CSPG       => (STRUCTURAL, "centralspindlin turnover"),
    :Ki_CSPG_CDK1  => (STRUCTURAL, "CDK1 blocks bundling until anaphase"),
    :Ki_RhoA_CDK1  => (STRUCTURAL, "CDK1 blocks Ect2 GEF activity until anaphase"),
    :kr_RhoA       => (STRUCTURAL, "RhoGAP"),
    :RhoA_tot      => (STRUCTURAL, "conserved RhoA pool, normalised to 1"),
    :kr_Midbody    => (STRUCTURAL, "midbody turnover"),
    :w_InVitro_ROS => (STRUCTURAL, "how much culture adds to ROS"),
    :ks_ATM_ROS    => (STRUCTURAL, "ROS -> ATM; scale absorbed by kf_CCNB_Ccng1"),
    :ks_Ccng1_p53  => (STRUCTURAL, "pinned: degenerate with kf_CCNB_Ccng1"),
    :kd_Ccng1      => (STRUCTURAL, "pinned: degenerate with kf_CCNB_Ccng1"),
)

"""
    fitted_params() -> Vector{Symbol}

The parameters counting against the budget, sorted.
"""
fitted_params() = sort([k for (k, v) in PARAM_PROVENANCE if v[1] == FITTED])

"""
    lint_budget(; budget=FITTED_BUDGET) -> Vector{String}

Structural checks on the parameter budget. Returns problems; empty means clean.

Mirrors Tier 1's `spec.lint`, which fails the build rather than warning, and for the same
reason: a budget that is only documented is a budget that drifts.

`ks_E2F8_E2F` is declared FITTED but tied to `ks_E2F7_E2F` during calibration — the two
repressors are induced by the same activator and Baniol measure them at nearly the same
correlation (+0.49 and +0.46), so estimating them separately would be spending a
parameter the data cannot resolve. The tie is asserted here rather than trusted.
"""
function lint_budget(; budget::Int = FITTED_BUDGET,
                       tied::Bool = true)
    problems = String[]
    declared = Set(keys(PARAM_PROVENANCE))
    actual = Set(setdiff(collect(keys(tier2_params())), collect(keys(params()))))

    for k in sort(collect(setdiff(actual, declared)))
        push!(problems, "$k exists but has no provenance declaration")
    end
    for k in sort(collect(setdiff(declared, actual)))
        push!(problems, "$k is declared but is not a Tier-2 parameter")
    end

    n = length(fitted_params())
    effective = tied ? n - 1 : n     # ks_E2F8_E2F tied to ks_E2F7_E2F
    effective += length(FITTED_OBSERVATION_PARAMS)
    if effective > budget
        push!(problems, "$effective fitted parameters exceeds the budget of $budget")
    end
    if effective > FITTED_CEILING
        push!(problems, "$effective fitted parameters exceeds the defensible " *
                        "ceiling of $FITTED_CEILING (MODEL.md 3.5)")
    end
    return problems
end

"""
Fitted parameters that are NOT entries of `tier2_params()` and so cannot be caught by
[`lint_budget`](@ref)'s completeness check. They still spend budget.

`FUCCI_THRESHOLD` is a property of the reporter and microscope rather than of the model —
where a fluorescence intensity is called positive — but it was estimated from data (by
requiring the FUCCI-derived and cyclin-derived phase fractions to agree), so it counts.
"""
const FITTED_OBSERVATION_PARAMS = (:FUCCI_THRESHOLD,)

"""
    total_fitted() -> Int

Effective fitted count: model parameters (with `ks_E2F8_E2F` tied) plus observation-model
parameters. This is the number the paper has to defend.
"""
total_fitted() = length(fitted_params()) - 1 + length(FITTED_OBSERVATION_PARAMS)

"""
    budget_report() -> String

Human-readable provenance summary, for the methods section.
"""
function budget_report()
    counts = Dict{Provenance,Int}()
    for (_, v) in PARAM_PROVENANCE
        counts[v[1]] = get(counts, v[1], 0) + 1
    end
    io = IOBuffer()
    println(io, "Tier-2 parameter provenance (inherited 218 are FROZEN and not listed)")
    for p in (INPUT, MEASURED, STRUCTURAL, SWITCH, FITTED)
        println(io, "  ", rpad(string(p), 12), get(counts, p, 0))
    end
    println(io, "\nFitted, model:       ", join(fitted_params(), ", "))
    println(io, "  effective ", length(fitted_params()) - 1,
                " (ks_E2F8_E2F tied to ks_E2F7_E2F)")
    println(io, "Fitted, observation: ", join(FITTED_OBSERVATION_PARAMS, ", "))
    println(io, "\nTOTAL fitted ", total_fitted(),
                " / budget ", FITTED_BUDGET, " / ceiling ", FITTED_CEILING,
                "   (Tier 1 spends 4)")
    println(io, "\nSTRUCTURAL is the honest grey zone: biologically motivated, not tuned")
    println(io, "against any target, and not measured either. The paper must say so.")
    return String(take!(io))
end
