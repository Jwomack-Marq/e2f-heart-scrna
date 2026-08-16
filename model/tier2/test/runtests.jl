using Test
using CmTier2
using Statistics
using Random

# ---------------------------------------------------------------------------
# Phase 0 gate: the ported model reproduces the published one.
#
# The inherited core in src/inherited/ is a byte-identical copy of the published
# `Cell_Cycle_Model/Current Iteration/model_files/`. Everything here pins behaviour so
# that when the Tier-2 cardiomyocyte modules land, "it still reduces to the published
# model" is a test result and not an assertion in a commit message.
#
# The source repo has no test suite at all, so these are also the first tests the
# inherited model has ever had.
# ---------------------------------------------------------------------------

const TOL = 1e-3   # hours, on period measurements

@testset "CmTier2 — Phase 0" begin

@testset "structure" begin
    @test length(state_names()) == 63
    @test length(state()) == 63
    @test length(params()) == 218
    @test allunique(state_names())

    # The state vector's component order must match state_names(), because diff_eqns.jl
    # destructures `u` positionally. A mismatch silently computes the wrong model.
    @test collect(keys(state())) == Symbol.(state_names())

    @test species_index("CCNB_CDK1") == 21
    @test species_index("CDT1") == 61
    @test_throws ErrorException species_index("NotASpecies")

    # The FUCCI observables TODO.md item 3 calls "the highest-value single feature".
    for s in ("CDT1", "Geminin", "Geminin_CDT1")
        @test s in state_names()
    end
    # Present, and the reason this model beats gg2009 as a Tier-2 base.
    for s in ("LMNA", "LMNAp", "PTTG1", "ppRB", "p21")
        @test s in state_names()
    end
end

@testset "baseline solves and oscillates" begin
    sol = solve_baseline(alpha = PUBLISHED_ALPHA)
    @test sol.retcode == :Success || string(sol.retcode) == "Success"

    idx = species_index("CCNB_CDK1")
    y = sol[idx, :]
    @test all(isfinite, y)
    @test minimum(y) >= -1e-8          # concentrations stay non-negative
    @test maximum(y) > 0.1             # the oscillation has real amplitude

    # A sustained limit cycle, not a decaying transient: amplitude late in the run
    # must match amplitude earlier. This is one of the validity criteria the source
    # repo's sensitivity pipeline applies but never asserts.
    early = [y[i] for i in eachindex(sol.t) if 1500 <= sol.t[i] <= 1700]
    late  = [y[i] for i in eachindex(sol.t) if 2300 <= sol.t[i] <= 2500]
    @test maximum(late) ≈ maximum(early) rtol = 0.02
end

@testset "published doubling times reproduce exactly" begin
    # α = 1.447 is what the publication figure scripts hardcode, and it reproduces the
    # committed experimental_vs_simulation_summary.csv (DMSO 28.1 h).
    @test peak_period(solve_baseline(alpha = 1.447)) ≈ 28.113 atol = TOL

    # α = 1.725 is implied by the default target in estimate_alpha_for_target_dt.jl.
    @test peak_period(solve_baseline(alpha = 1.725)) ≈ 23.582 atol = TOL

    # α = 1.60 is what the varyingVolo figure variant uses.
    @test peak_period(solve_baseline(alpha = 1.60)) ≈ 25.425 atol = TOL

    # α = 2.3 is the value sitting in parameters.jl. It corresponds to no published
    # figure. Pinned deliberately: if someone reconciles the four α values, this test
    # fails and forces the decision to be recorded rather than absorbed.
    @test peak_period(solve_baseline(alpha = PARAMFILE_ALPHA)) ≈ 17.685 atol = TOL
end

@testset "period scales as 1/alpha" begin
    # α multiplies the whole RHS, so this must hold exactly. It is what makes the four
    # competing α values a pure rescaling rather than four different models.
    base = peak_period(solve_baseline(alpha = 1.0))
    for a in (1.447, 1.725, 2.3)
        @test peak_period(solve_baseline(alpha = a)) ≈ base / a rtol = 1e-3
    end
end

@testset "the two period measurements agree" begin
    # doubling_time() is the threshold-crossing method from the α helper;
    # peak_period() is the peak-detection method from the figure scripts. Independent
    # implementations, so agreement is a real check on both.
    sol = solve_baseline(alpha = PUBLISHED_ALPHA)
    @test doubling_time(sol) ≈ peak_period(sol) atol = 5e-3
end

@testset "drug dose-response matches the committed CSV" begin
    # Targets read from Cell_Cycle_Model/Current Iteration/experimental_vs_simulation_summary.csv,
    # measured post-dose in (2000, 2500) exactly as publication_fig_12 does.
    win = (2000.0, 2500.0)
    cases = [
        # (vola nM, abe µM, published Sim_Doubling_Time_hr)
        (1.0, 0.0,   28.1),
        (1.0, 0.001, 28.1),
        (1.0, 0.01,  28.2),
        (1.0, 0.1,   40.7),
    ]
    for (vola, abe, published) in cases
        sol = solve_drug(vola_nM = vola, abe_uM = abe, alpha = PUBLISHED_ALPHA)
        got = peak_period(sol; window = win)
        @test got ≈ published atol = 0.1
    end

    # The dose-response has the right sign and the 0.1 µM step is the big one.
    dt_low  = peak_period(solve_drug(vola_nM = 1.0, abe_uM = 0.001); window = win)
    dt_high = peak_period(solve_drug(vola_nM = 1.0, abe_uM = 0.1);   window = win)
    @test dt_high > dt_low
end

@testset "known equation defects are pinned, not silently inherited" begin
    # These are recorded so Phase 2 cannot fix them by accident and call it a wash.
    # See the plan's "What exactly is missing, verified against the source".
    src = read(joinpath(@__DIR__, "..", "src", "inherited", "diff_eqns.jl"), String)

    # 1. CDT1/Geminin synthesis is constitutive despite the `_E2F` parameter names, so
    #    the FUCCI observables are currently decoupled from E2F. Phase 2 connects them.
    @test occursin(r"d\.CDT1\s*=\s*\(\s*\n?\s*ks_CDT1_E2F\s*\n", src) ||
          occursin("ks_CDT1_E2F\n", src)
    @test !occursin("ks_CDT1_E2F*E2F", src)
    @test !occursin("ks_Geminin_E2F*E2F", src)

    # 2. d.p53p is missing the *p53 factor that d.p53's matching loss term carries, so
    #    p53 <-> p53p is not mass-balanced. Harmless while the DDR arm is dormant
    #    (kf_ATMp = 0), which this asserts, but Phase 2 turns that arm on.
    @test occursin("d.p53p = ((kf_p53p*Chk2p)-kr_p53p*PPase*p53p) * α", src)
    @test params().kf_ATMp == 0.0
end

@testset "E2F has exactly six downstream consumers" begin
    # The E2F sub-family split in Phase 2 must rewire each of these and nothing else.
    # If the inherited model gains an E2F consumer, the split silently misses it.
    src = read(joinpath(@__DIR__, "..", "src", "inherited", "diff_eqns.jl"), String)
    consumers = ["ks_CCNE_E2F*E2F", "ks_CDC25A_E2F*E2F", "ks_CCNA_E2F*E2F",
                 "ks_PLK1 * E2F", "ks_PTTG1 * E2F", "ks_EMI1_E2F*E2F"]
    for c in consumers
        @test occursin(c, src)
    end
end

@testset "the cardiomyocyte layer is genuinely absent" begin
    # Phase 2 adds these. Asserting their absence now makes the diff auditable and
    # stops a half-ported module from going unnoticed.
    names = state_names()
    for absent in ("Ect2", "RhoA", "Anillin", "AurKB", "Midbody", "p27", "Ccng1")
        @test !(absent in names)
    end
end

end # testset


# ---------------------------------------------------------------------------
# Phase 1: observables, discrete landmarks, fate classification.
#
# Adds no biology and no parameters — it only measures the inherited model. The gate is
# whether the measured durations are commensurate with the five `duration` rows in
# ../cmcycle/data/cmcycle_targets.csv, which Tier 1 structurally cannot produce.
# ---------------------------------------------------------------------------

"""Settled-limit-cycle window: ~14 cycles at the published 28.11 h period."""
const W = (1800.0, 2200.0)

@testset "CmTier2 — Phase 1" begin

@testset "events fire once per cycle, in biological order" begin
    sol, log = solve_with_events(alpha = PUBLISHED_ALPHA)
    l = trim(log, W)

    # One of each landmark per cycle. The window spans ~14 cycles at 28.11 h.
    n = length(l.restriction)
    @test n >= 10
    for v in (l.s_entry, l.neb, l.anaphase, l.mitotic_exit, l.geminin_on)
        @test abs(length(v) - n) <= 1        # allow one boundary-clipped event
    end

    # Within a cycle: commit -> S entry -> NEB -> anaphase -> exit.
    cycles = classify_cycles(log; window = W)
    @test length(cycles) >= 10
    for c in cycles
        @test c.neb !== nothing && c.anaphase !== nothing && c.exit !== nothing
        @test c.s_entry < c.neb < c.exit
        @test c.anaphase <= c.exit           # securin destroyed before MPF is gone
        @test c.mitosis > 0
    end

    # Deterministic limit cycle: every traversal identical.
    @test all(c -> c.fate === :MitoticCompletion, cycles)
    sg2m = [c.sg2m_cyclin for c in cycles]
    @test maximum(sg2m) - minimum(sg2m) < 1e-3
end

@testset "recording events does not perturb the trajectory" begin
    # The affects push onto a log and never touch u or p, so the solution must be the
    # same ODE solution. If this fails, every number measured here is suspect.
    #
    # It cannot be asserted at default tolerances: a ContinuousCallback forces the
    # solver to stop at each event root, which changes step selection, so the two solves
    # differ by ~2e-3 relative -- the solver's own reltol, not a perturbation. The
    # discriminating test is that the difference CONVERGES TO ZERO with tolerance, which
    # a genuine perturbation would not. Measured: 2e-3 (default) -> 5.9e-5 (1e-6) ->
    # 3.3e-9 (1e-10).
    diffs = Float64[]
    for tol in (1e-6, 1e-10)
        plain = solve_baseline(alpha = PUBLISHED_ALPHA, reltol = tol, abstol = tol)
        withcb, _ = solve_with_events(alpha = PUBLISHED_ALPHA, reltol = tol, abstol = tol)
        worst = 0.0
        for s in ("CCNB_CDK1", "CCNE_CDK2", "CDT1", "Geminin", "LMNAp")
            i = species_index(s)
            for t in (2000.0, 2137.5)
                a, b = withcb(t)[i], plain(t)[i]
                worst = max(worst, abs(a - b) / max(abs(b), 1e-12))
            end
        end
        push!(diffs, worst)
        @test peak_period(withcb) ≈ peak_period(plain) atol = 1e-3
    end
    @test diffs[1] < 1e-3          # already tight at reltol 1e-6
    @test diffs[2] < 1e-7          # and vanishes at 1e-10
    @test diffs[2] < diffs[1]      # converging, i.e. discretisation not perturbation
end

@testset "fate classification is total and exclusive" begin
    _, log = solve_with_events(alpha = PUBLISHED_ALPHA)
    cycles = classify_cycles(log; window = W)
    for c in cycles
        @test c.fate in FATES_PHASE1
    end
    s = fate_summary(log; window = W)
    @test sum(values(s[:counts])) == s[:n_cycles]
    @test sum(values(s[:fractions])) ≈ 1.0
end

@testset "quiescence is detected, and is dose-monotone" begin
    # CDK4/6 inhibition blocks Rb phosphorylation, so the restriction point is never
    # passed. This is the negative control for the classifier: a fate model that cannot
    # produce a quiescent cell cannot produce Tier 1's largest fate (90.35 %).
    w = (800.0, 1200.0)
    counts = Int[]
    for abe in (0.0, 0.1, 1.0, 10.0)
        _, log = solve_with_events(alpha = PUBLISHED_ALPHA, tspan = (0.0, 1200.0),
                                   con_ABE = abe)
        push!(counts, length(classify_cycles(log; window = w)))
        if abe >= 1.0
            @test quiescent(log; window = w)
        else
            @test !quiescent(log; window = w)
        end
    end
    @test issorted(counts, rev = true)   # more drug, never more cycles
    @test counts[end] == 0

    # A quiescent trajectory yields no cycles and the untouched initial bookkeeping.
    _, log = solve_with_events(alpha = PUBLISHED_ALPHA, tspan = (0.0, 1200.0),
                               con_ABE = 10.0)
    s = fate_summary(log; window = w)
    @test s[:quiescent]
    @test s[:n_cycles] == 0
    @test s[:book].cells == 1.0
end

@testset "bookkeeping counts cells, nuclei and ploidy" begin
    _, log = solve_with_events(alpha = PUBLISHED_ALPHA)
    cycles = classify_cycles(log; window = W)

    # n rounds of division from one cell.
    b = bookkeep(cycles; assume_abscission = true)
    @test b.cells ≈ 2.0^length(cycles)
    @test b.nuclei == 1.0
    @test b.dna_content == 2.0            # each daughter is diploid again

    # Same cycles with abscission failing: one cell, one nucleus added per round.
    bb = bookkeep(cycles; assume_abscission = false)
    @test bb.cells == 1.0
    @test bb.nuclei ≈ 1.0 + length(cycles)

    # Polyploidization keeps doubling DNA in a single nucleus.
    poly = [Cycle(:Polyploidization, 0.0, nothing, nothing, nothing, nothing,
                  nothing, nothing, nothing, nothing) for _ in 1:3]
    bp = bookkeep(poly)
    @test bp.cells == 1.0
    @test bp.nuclei == 1.0
    @test bp.dna_content ≈ 16.0           # 2C -> 4C -> 8C -> 16C
end

@testset "phase durations resolve and are self-consistent" begin
    sol = solve_baseline(alpha = PUBLISHED_ALPHA)
    pt = phase_times(sol)
    @test pt !== nothing
    @test pt.g1 + pt.s + pt.g2 + pt.m ≈ pt.cell_cycle rtol = 1e-6
    @test all(>(0), (pt.g1, pt.s, pt.g2, pt.m))
    @test 0 < pt.g1_percent < 100

    # phase_times reads boundaries off peak INDICES, so its cycle length is quantised to
    # the save grid, and CCNE peaks in a slowly-varying stretch where the adaptive solver
    # takes long steps. On the default grid it disagrees with peak_period by 0.097 h
    # (0.34 %). That is quantisation, not disagreement about the period: refining the
    # grid drives it to zero -- 0.097 h (auto) -> 0.033 h (saveat 0.05) -> 0.0067 h
    # (saveat 0.01). Asserted as convergence so the two measurements cannot silently
    # drift apart for a real reason.
    gaps = [abs(phase_times(s).cell_cycle - peak_period(s)) for s in
            (sol,
             solve_baseline(alpha = PUBLISHED_ALPHA, saveat = 0.05),
             solve_baseline(alpha = PUBLISHED_ALPHA, saveat = 0.01))]
    @test gaps[1] < 0.15
    @test gaps[3] < 0.02
    @test gaps[3] < gaps[1]

    # Absolute peak criteria cannot serve all four species: free CCNE peaks at ~0.13, so
    # the height = 0.2 used for CCNB_CDK1 (peak 0.66) finds no CCNE peak at all and the
    # G1/S boundary is never located. Pinned so the relative default is not "simplified"
    # back to an absolute one.
    @test maximum(sol[species_index("CCNE"), :]) < 0.2
    @test maximum(sol[species_index("CCNB_CDK1"), :]) > 0.5
end

@testset "durations are commensurate with the published measurements" begin
    _, log = solve_with_events(alpha = PUBLISHED_ALPHA)
    s = fate_summary(log; window = W)

    # Murganti Fig 2E S/G2/M: 16.38 (division), 17.29 (binucleation), 24.50 h
    # (polyploidization). Baniol P0 mouse 15.1 +/- 4.0 SD. Nothing here is fitted -- the
    # inherited model is a generic proliferating line -- so the test is order of
    # magnitude and range, not agreement.
    @test 10.0 < s[:mean_sg2m_cyclin] < 30.0

    # Mitosis. Tier 1's preflight derives ~4.2 h for cardiomyocytes; a generic cycling
    # line is faster. Bounded well away from both zero and the whole cycle.
    @test 0.5 < s[:mean_mitosis] < 6.0
    @test s[:mean_mitosis] < 0.2 * peak_period(solve_baseline(alpha = PUBLISHED_ALPHA))
end

@testset "RETRACTED: the FUCCI 'defect' was the 0.05 plotting default" begin
    # This testset previously asserted a KNOWN FAILURE -- that the inherited FUCCI layer
    # had no G1/S double-positive state. Phase 2 showed that was an artefact of the
    # cutoff, not a property of the model, so the finding is retracted and the tests are
    # inverted to pin the correct behaviour.
    #
    # The cutoff is a reporter/microscope property, not a model property. Calibrating it
    # by requiring the FUCCI-derived and cyclin-derived phase fractions to agree -- two
    # measurements sharing no equations -- gives 0.02 with a sharp optimum, and there the
    # licensing layer is excellent.
    sol = solve_baseline(alpha = PUBLISHED_ALPHA)
    pt = phase_times(sol)
    g1 = pt.g1 / pt.cell_cycle
    sg2m = (pt.s + pt.g2 + pt.m) / pt.cell_cycle

    f = fucci_fractions(sol; window = W)          # calibrated cutoff
    @test f[:G0G1] + f[:G1S] ≈ g1   atol = 0.02   # Cdt1-positive == G1
    @test f[:SG2M]           ≈ sg2m atol = 0.02   # geminin-positive == S+G2+M
    @test f[:G1S] > 0.05                          # a real double-positive exists
    @test f[:negative] < 0.01                     # and essentially no double-negatives
    @test sum(values(f)) ≈ 1.0

    # At the old plotting default the same model looks broken. Pinned so nobody
    # "restores" it and re-derives the retracted conclusion.
    fp = fucci_fractions(sol; window = W, threshold = PUBLISHED_FUCCI_THRESHOLD)
    @test fp[:G1S] == 0.0
    @test fp[:negative] > 0.30

    # The optimum is sharp, i.e. the cutoff is genuinely identified rather than a range.
    err(thr) = (x = fucci_fractions(sol; window = W, threshold = thr);
                abs(x[:G0G1] + x[:G1S] - g1) + abs(x[:SG2M] - sg2m) + x[:negative])
    @test err(FUCCI_THRESHOLD) < 0.01
    @test err(0.025) > 5 * err(FUCCI_THRESHOLD)
    @test err(0.015) > 5 * err(FUCCI_THRESHOLD)
end

@testset "FUCCI-timed duration at the calibrated cutoff" begin
    # With the cutoff calibrated, the FUCCI-timed S/G2/M is 12.80 h -- against Murganti
    # Fig 2E's 16.38 h and Baniol's 15.1 +/- 4.0 SD for P0 mouse, i.e. inside one SD of
    # the mouse measurement with nothing fitted to either. At the old 0.05 default the
    # same quantity read 2.06 h.
    #
    # It should sit between the phase_times S+G2+M (10.73 h, cyclin peaks) and the
    # CycE-crossing measure (18.95 h): geminin appears before CycA peaks but after CycE
    # crosses, so the three bracket each other in a fixed order.
    _, log = solve_with_events(alpha = PUBLISHED_ALPHA)
    s = fate_summary(log; window = W)
    sol = solve_baseline(alpha = PUBLISHED_ALPHA)
    pt = phase_times(sol)

    @test 10.0 < s[:mean_sg2m_fucci] < 20.0
    @test s[:mean_sg2m_fucci] > pt.s + pt.g2 + pt.m
    @test s[:mean_sg2m_fucci] < s[:mean_sg2m_cyclin]
    @test abs(s[:mean_sg2m_fucci] - 15.1) < 4.0      # within Baniol's SD
end

@testset "fucci_state truth table" begin
    @test fucci_state(0.9, 0.9) === :G1S
    @test fucci_state(0.9, 0.0) === :G0G1
    @test fucci_state(0.0, 0.9) === :SG2M
    @test fucci_state(0.0, 0.0) === :negative
    @test fucci_state(FUCCI_THRESHOLD, FUCCI_THRESHOLD) === :negative  # strict >
end

@testset "total_pools restores GemininT" begin
    # post_process_cc in the source repo has the GemininT block commented out while
    # plot_fucci_backgrounds reads solution_df.GemininT, so that path would raise.
    sol = solve_baseline(alpha = PUBLISHED_ALPHA)
    p = total_pools(sol)
    @test all(p.GemininT .≈ sol[species_index("Geminin"), :] .+
                            sol[species_index("Geminin_CDT1"), :])
    @test maximum(p.GemininT) > maximum(sol[species_index("Geminin"), :])
    for k in (:WEE1T, :PLK1T, :APCCT)
        @test all(isfinite, getfield(p, k))
    end
end

end # testset



# ---------------------------------------------------------------------------
# Phase 2: the cardiomyocyte modules.
#
# Step 1 (coupling Cdt1/geminin synthesis to E2F) was TESTED AND REJECTED -- see below.
# Step 2 (the E2F sub-family split) is wired and inert by default.
# ---------------------------------------------------------------------------

@testset "CmTier2 — Phase 2" begin

@testset "THE GATE: disabled Tier 2 is bit-exactly the published model" begin
    # The strongest form of the reduction property. tier2DiffEq! CALLS the inherited RHS
    # rather than copying it, and every Tier-2 term carries a factor that is zero at
    # default parameters -- so the reduction is structural, not a tolerance.
    #
    # Checked on random states rather than on a trajectory, because a trajectory only
    # visits the limit cycle and would miss a correction that is nonzero elsewhere.
    rng_state = 0
    worst = 0.0
    for trial in 1:100
        u63, u66 = state(), tier2_state()
        for (k, s) in enumerate(state_names())
            rng_state = (1103515245 * rng_state + 12345) % 2147483648
            v = 0.8 * rng_state / 2147483648
            u63[Symbol(s)] = v
            u66[Symbol(s)] = v
        end
        d63, d66 = similar(u63), similar(u66)
        modelDiffEq!(d63, u63, params(), 0.0)
        tier2DiffEq!(d66, u66, tier2_params(), 0.0)
        for i in 1:63
            worst = max(worst, abs(d66[i] - d63[i]))
        end
    end
    @test worst == 0.0          # bit-exact, not approximately

    # And the same at the trajectory level, where solver step selection differs slightly
    # because the error norm now spans 66 components rather than 63.
    plain = solve_baseline(alpha = PUBLISHED_ALPHA, reltol = 1e-10, abstol = 1e-10)
    t2    = solve_tier2(alpha = PUBLISHED_ALPHA, reltol = 1e-10, abstol = 1e-10)
    for s in ("CCNB_CDK1", "CCNE_CDK2", "E2F", "CDT1", "Geminin")
        i = species_index(s)
        @test t2(2000.0)[i] ≈ plain(2000.0)[i] rtol = 1e-5
    end
    @test peak_period(t2) ≈ peak_period(plain) atol = 1e-3
    for s in TIER2_SPECIES
        @test abs(t2(2000.0)[species_index(s)]) < 1e-30
    end
end

@testset "extended vectors keep the inherited layout" begin
    @test length(tier2_state()) == 63 + length(TIER2_SPECIES)
    # Tripwire: bumps by design as each step lands. Currently 41 Tier-2 params.
    @test length(tier2_params()) == 218 + 41
    @test tier2_state_names()[1:63] == state_names()
    # New components MUST come last: diff_eqns.jl destructures positionally.
    @test tier2_state_names()[64:end] == collect(TIER2_SPECIES)
    @test length(TIER2_SPECIES) == 10      # tripwire: bumps by design as steps land
    @test collect(keys(tier2_params()))[1:218] == collect(keys(params()))
    for s in TIER2_SPECIES
        @test species_index(s) > 63
    end
    # Every enable parameter defaults to zero -- that IS the reduction property.
    for k in TIER2_ENABLE_PARAMS
        @test getproperty(tier2_params(), k) == 0.0
    end
end

@testset "e2f_repression is well behaved" begin
    @test e2f_repression(0, 0, 0, 0.1, 0.2) == 1.0        # unrepressed
    @test e2f_repression(0, 1, 0, 0.1, 0.2) < 1.0
    @test e2f_repression(0, 0, 1, 0.1, 0.2) < 1.0
    @test e2f_repression(1, 0, 0, 0.1, 0.2) < 1.0
    # E2F7 and E2F8 act as one pool: only their sum matters.
    @test e2f_repression(0, 0.3, 0.1, 0.1, 0.2) == e2f_repression(0, 0.1, 0.3, 0.1, 0.2)
    # E2F6 represses more weakly at equal concentration (larger Ki).
    @test e2f_repression(0.5, 0, 0, 0.1, 0.2) > e2f_repression(0, 0.5, 0, 0.1, 0.2)
    # Bounded in (0, 1] and monotone.
    @test 0 < e2f_repression(10, 10, 10, 0.1, 0.2) <= 1
    @test e2f_repression(0, 2, 0, 0.1, 0.2) < e2f_repression(0, 1, 0, 0.1, 0.2)
end

@testset "the E2F split reproduces Baniol's signs" begin
    # Baniol measured E2f1~E2f7 = +0.49, E2f1~E2f8 = +0.46, E2f8~Ccna2 = -0.56.
    #
    # Scored as SIGNS, not magnitudes, and deliberately so: those are Spearman
    # correlations across single cells at Smart-seq2 depth, where dropout attenuates
    # |r| substantially. A deterministic trajectory has no measurement error, so it will
    # always give larger |r|. Comparing magnitudes without an attenuation correction
    # would be scoring the noise model, not the biology.
    en = (ks_E2F7_E2F = 0.20, ks_E2F8_E2F = 0.20, Ki_E2F78 = 5.0)
    sol = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en, saveat = 0.25)
    t = sol.t
    w = findall(x -> 1900 <= x <= 2180, t)
    g(s) = sol[species_index(s), :][w]
    function pearson(a, b)
        am, bm = sum(a)/length(a), sum(b)/length(b)
        den = sqrt(sum((a .- am).^2) * sum((b .- bm).^2))
        den < 1e-14 ? NaN : sum((a .- am) .* (b .- bm)) / den
    end

    @test pearson(g("E2F"), g("E2F7")) > 0      # E2F1 induces E2F7
    @test pearson(g("E2F"), g("E2F8")) > 0      # and E2F8
    @test pearson(g("E2F8"), g("CCNA")) < 0     # E2F8 is G1/S-restricted
    @test maximum(g("E2F7")) > 0                # the repressors are actually expressed
    @test maximum(g("E2F8")) > 0
    # E2F8 turns over faster than E2F7 (CycA clearance), so it accumulates less.
    @test maximum(g("E2F8")) < maximum(g("E2F7"))
end

@testset "in-silico E2f7/E2f8 knockdown shortens the cycle" begin
    # The direction the lab's own data shows: cycling cardiomyocytes at P7 are 31.6 % in
    # the KO against 25.6 % WT, i.e. losing the repressors means MORE cycling. Here that
    # appears as a shorter period. Nothing was fitted to the KO data.
    base = (ks_E2F7_E2F = 0.20, ks_E2F8_E2F = 0.20, Ki_E2F78 = 5.0)
    per(over) = peak_period(solve_tier2(alpha = PUBLISHED_ALPHA,
                                        enable = merge(base, over),
                                        tspan = (0.0, 3000.0));
                            window = (2400.0, 2900.0))
    wt     = per(NamedTuple())
    kd7    = per((ks_E2F7_E2F = 0.0,))
    kd8    = per((ks_E2F8_E2F = 0.0,))
    double = per((ks_E2F7_E2F = 0.0, ks_E2F8_E2F = 0.0))

    @test all(isfinite, (wt, kd7, kd8, double))
    @test kd7 < wt && kd8 < wt                 # each knockdown speeds the cycle
    @test double < kd7 && double < kd8         # the double goes furthest
    # Removing both repressors must land exactly on the published model.
    @test double ≈ peak_period(solve_baseline(alpha = PUBLISHED_ALPHA)) atol = 0.05

    # Real epistasis: the double is not the sum of the singles. Tier 1 found the same
    # qualitative result on Ect2, which matters because the lab's data IS a double KO.
    @test abs((wt - double) - ((wt - kd7) + (wt - kd8))) > 0.5
end

@testset "REJECTED: coupling Cdt1/geminin synthesis to E2F" begin
    # ks_CDT1_E2F and ks_Geminin_E2F are constants despite their names -- a real naming
    # defect, pinned in Phase 0. Connecting them to E2F was the planned step 1. It was
    # implemented, measured, and REJECTED: it makes the FUCCI structure worse, not
    # better. Recorded here rather than deleted, following the project's practice of
    # committing negative results.
    #
    # Scored by how well the FUCCI-derived phase fractions match the cyclin-derived
    # ones, each variant at ITS OWN best cutoff so the comparison is fair:
    #
    #   published (off)      err 0.0006
    #   gem 1.0              err 0.4663
    #   cdt1 1.0             err 0.0973
    #   gem + cdt1 1.0       err 0.1524
    #   gem + cdt1 0.5       err 0.0402
    #
    # Mechanism of the failure: geminin's peak moves from +2.0 h (mitosis, correct for
    # mAG) to +16.5 h (late G1) and its amplitude falls from 0.133 to 0.043, because E2F
    # is a sharp late-G1 spike while geminin needs to accumulate across S/G2/M.
    #
    # kd_CDT1_CDK2 was also tried, to clear Cdt1 at S onset -- CDK2 phosphorylation of
    # Cdt1 licensing SCF-Skp2 is real biology, and Cdt1's only route here is SCF, which
    # peaks at mitosis rather than at G1/S. It is too blunt: CCNE_CDK2 + CCNA_CDK2 never
    # falls below ~0.07, so it acts as near-constant degradation and erases Cdt1
    # entirely (all four FUCCI fractions collapse to double-negative by kd = 5). Doing
    # it properly needs an S-phase marker the inherited model does not have -- there is
    # no DNA replication variable. Left wired, defaulted off, for Phase 3.
    W2 = (1800.0, 2200.0)
    # Fixed save grid: phase_times reads boundaries off peak INDICES, so on the adaptive
    # grid this score wanders with step selection (0.0252 vs 0.0068). The variants are
    # unaffected, so the verdict never depended on it -- but the baseline number did.
    function fit_err(en)
        sol = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en, saveat = 0.05)
        pt = phase_times(sol)
        pt === nothing && return Inf
        g1 = pt.g1 / pt.cell_cycle
        sg = (pt.s + pt.g2 + pt.m) / pt.cell_cycle
        best = Inf
        for thr in 0.004:0.002:0.05
            f = fucci_fractions(sol; window = W2, threshold = thr)
            best = min(best, abs(f[:G0G1] + f[:G1S] - g1) + abs(f[:SG2M] - sg) + f[:negative])
        end
        return best
    end
    off = fit_err(NamedTuple())
    @test off < 0.015                                       # published is near-perfect
    @test fit_err((w_Geminin_E2F = 1.0,)) > 10 * off        # and every variant is worse
    @test fit_err((w_CDT1_E2F = 1.0,)) > off
    @test fit_err((w_Geminin_E2F = 0.5, w_CDT1_E2F = 0.5)) > off

    # The blunt-clearance failure, pinned. The published model has essentially no
    # double-negative population at the calibrated cutoff; adding CDK2-driven Cdt1
    # clearance creates a large one, because the clearance is on whenever CDK2 is, which
    # is most of the cycle.
    base_neg = fucci_fractions(solve_tier2(alpha = PUBLISHED_ALPHA); window = W2)[:negative]
    f = fucci_fractions(solve_tier2(alpha = PUBLISHED_ALPHA,
                                    enable = (w_Geminin_E2F = 1.0, w_CDT1_E2F = 1.0,
                                              kd_CDT1_CDK2 = 5.0)); window = W2)
    @test base_neg < 0.01
    @test f[:negative] > 0.4
    @test f[:negative] > 50 * max(base_neg, 1e-4)
end

end # testset


# ---------------------------------------------------------------------------
# Phase 2 step 3: the cytokinesis arm.
# ---------------------------------------------------------------------------

const E2F_ON = (ks_E2F7_E2F = 0.20, ks_E2F8_E2F = 0.20, Ki_E2F78 = 5.0)
const W3 = (1900.0, 2180.0)

@testset "CmTier2 — Phase 2 step 3 (cytokinesis)" begin

@testset "the arm is inert until switched on" begin
    # The reduction property must survive every new species. The total count is a
    # tripwire in the layout testset; here we only care that the arm's six are present.
    for s in ("Ect2", "RhoA", "Centralspindlin", "AurKB", "Anillin", "Midbody")
        @test s in TIER2_SPECIES
        @test species_index(s) > 63
    end
    sol = solve_tier2(alpha = PUBLISHED_ALPHA)
    for s in ("Ect2", "RhoA", "Midbody", "Anillin", "Centralspindlin", "AurKB")
        @test abs(sol(2000.0)[species_index(s)]) < 1e-30
    end
end

@testset "the arm runs and produces a midbody" begin
    en = merge(E2F_ON, CYTOKINESIS_ON)
    sol, log = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en,
                           record_events = true, cytokinesis = true)
    for s in ("Ect2", "AurKB", "Centralspindlin", "RhoA", "Anillin", "Midbody")
        v = sol[species_index(s), :]
        @test all(isfinite, v)
        @test minimum(v) >= -1e-8
        @test maximum(v) > 0
    end
    l = trim(log, W3)
    @test length(l.abscission) > 0
    # One abscission per mitosis: the midbody is not re-crossing within a cycle.
    @test abs(length(l.abscission) - length(l.neb)) <= 1

    cycles = classify_cycles(log; window = W3, cytokinesis = true)
    @test !isempty(cycles)
    @test all(c -> c.fate === :Division, cycles)
    @test all(c -> c.abscission !== nothing, cycles)
end

@testset "Ect2 is rate-limiting for division, and only for division" begin
    # Tier 1's central claim, now mechanistic: knocking Ect2 down must remove division
    # while leaving S-phase entry and mitotic entry untouched. A model in which Ect2
    # knockdown also reduced entry would be describing a general toxicity, not a
    # cytokinesis-specific block.
    function probe(over)
        en = merge(merge(E2F_ON, CYTOKINESIS_ON), over)
        sol, log = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en,
                               record_events = true, cytokinesis = true)
        cyc = classify_cycles(log; window = W3, cytokinesis = true)
        l = trim(log, W3)
        w = findall(x -> W3[1] <= x <= W3[2], sol.t)
        (div = count(c -> c.fate === :Division, cyc),
         bi  = count(c -> c.fate === :Binucleation, cyc),
         s_entries = length(l.s_entry), mitoses = length(l.neb),
         midbody = maximum(sol[species_index("Midbody"), :][w]))
    end

    wt = probe(NamedTuple())
    kd = probe((ks_Ect2_E2F = 0.0,))

    @test wt.div > 0 && wt.bi == 0            # WT divides
    @test kd.div == 0                          # knockdown removes 100 % of division
    @test kd.bi > 0                            # and converts it to binucleation
    @test kd.s_entries == wt.s_entries         # entry untouched
    @test kd.mitoses    == wt.mitoses          # mitotic entry untouched
    @test kd.midbody < 1e-12

    # Graded, not just on/off: the midbody scales with Ect2. That is what will turn into
    # graded fate FRACTIONS once Phase 3 runs a heterogeneous population -- a single
    # deterministic cell can only ever give one fate.
    mids = [probe((ks_Ect2_E2F = 0.60 * f,)).midbody for f in (0.0, 0.1, 0.3, 0.5, 1.0)]
    @test issorted(mids)
    @test mids[end] > 5 * mids[2]
end

@testset "RhoA is obligatory by construction" begin
    # d.Midbody has exactly ONE production term with RhoA as a factor, so no parameter
    # setting can make a midbody without it. Tier 1 hit the opposite arrangement four
    # times; the worst OR'd a second route on and made it structurally impossible for
    # Ect2 to be rate-limiting -- the model could not express its own central claim.
    src = read(joinpath(@__DIR__, "..", "src", "tier2_model.jl"), String)
    code = filter(l -> !startswith(strip(l), "#"), split(src, '\n'))
    midbody_lines = filter(l -> occursin("d.Midbody", l), code)
    @test length(midbody_lines) == 1
    @test occursin("kf_Midbody * RhoA * Anillin", midbody_lines[1])

    # And behaviourally: every upstream component is individually necessary.
    function midbody_max(over)
        en = merge(merge(E2F_ON, CYTOKINESIS_ON), over)
        sol = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en)
        w = findall(x -> W3[1] <= x <= W3[2], sol.t)
        maximum(sol[species_index("Midbody"), :][w])
    end
    @test midbody_max(NamedTuple()) > 0.05
    for ko in ((kf_RhoA = 0.0,), (ks_Anln_E2F = 0.0,), (ks_AurKB_CDK1 = 0.0,),
               (ks_CSPG = 0.0,), (ks_Ect2_E2F = 0.0,))
        @test midbody_max(ko) < 1e-12
    end
end

@testset "the four-fate partition, and Phase 1 labelling is preserved" begin
    @test FATES == (:Quiescent, :Division, :Binucleation, :Polyploidization)

    # Without the arm, division and binucleation are indistinguishable and the honest
    # Phase 1 label must survive. Reporting :Binucleation just because no midbody
    # formed would be a guess dressed as a result.
    _, log = solve_with_events(alpha = PUBLISHED_ALPHA)
    plain = classify_cycles(log; window = W)
    @test all(c -> c.fate === :MitoticCompletion, plain)
    @test all(c -> c.abscission === nothing, plain)

    # With the arm, the same machinery resolves the split.
    en = merge(E2F_ON, CYTOKINESIS_ON)
    _, log2 = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en,
                          record_events = true, cytokinesis = true)
    split = classify_cycles(log2; window = W3, cytokinesis = true)
    @test all(c -> c.fate in FATES, split)
    @test !any(c -> c.fate === :MitoticCompletion, split)

    s = fate_summary(log2; window = W3, cytokinesis = true)
    @test sum(values(s[:counts])) == s[:n_cycles]
    @test sum(values(s[:fractions])) ≈ 1.0
end

@testset "bookkeeping follows the model's own fate call" begin
    # With the arm on, assume_abscission must not be consulted: the outcome is decided.
    divs = [Cycle(:Division, 0.0, nothing, nothing, nothing, nothing, 1.0,
                  nothing, nothing, nothing) for _ in 1:4]
    b1 = bookkeep(divs; assume_abscission = true)
    b2 = bookkeep(divs; assume_abscission = false)
    @test b1.cells == b2.cells == 16.0        # 2^4, either way
    @test b1.nuclei == 1.0
    @test b1.dna_content == 2.0

    # Binucleation: one cell, a nucleus per failed round.
    bis = [Cycle(:Binucleation, 0.0, nothing, nothing, nothing, nothing, nothing,
                 nothing, nothing, nothing) for _ in 1:3]
    bb = bookkeep(bis)
    @test bb.cells == 1.0
    @test bb.nuclei == 4.0                    # 1 + 3
    @test bb.dna_content == 2.0

    # A mixed history: 2 divisions then 1 binucleation.
    mixed = vcat(divs[1:2], bis[1:1])
    bm = bookkeep(mixed)
    @test bm.cells == 4.0
    @test bm.nuclei == 2.0
end

end # testset


# ---------------------------------------------------------------------------
# Phase 2 step 4: the maturation axis.
# ---------------------------------------------------------------------------

@testset "CmTier2 — Phase 2 step 4 (maturation)" begin

@testset "contexts come from Tier 1, not a copy" begin
    # If a context drifts between tiers that is a bug, not a configuration -- cross-tier
    # comparison is only meaningful if both tiers mean the same cell.
    @test isfile(CMFATE_MANIFEST)
    @test occursin("cmcycle", CMFATE_MANIFEST)
    c = contexts()
    @test length(c) == 6
    for name in ("hipsc_cm", "mncm_invitro", "mouse_p0_invivo", "mouse_p1_invivo",
                 "mouse_p7_invivo", "adult")
        @test haskey(c, name)
    end
    @test maturation("hipsc_cm") == 0.12
    @test maturation("mouse_p0_invivo") == 0.30
    @test maturation("mouse_p7_invivo") == 0.55
    @test maturation("adult") == 0.95
    @test context_names()[1] == "hipsc_cm"      # sorted by M
    @test context_names()[end] == "adult"
    @test_throws ErrorException maturation("not_a_context")
end

@testset "the couplings are measured, and cost one parameter not two" begin
    # M is z-scored within one 285-cell dataset and has no absolute cross-system scale,
    # so a correlation between two z-scored quantities fixes SIGNS and their RATIO but
    # not the absolute slopes. Hence one shared gain with the ratio welded in.
    @test MATURATION_SLOPE_ECT2 == -0.563     # Baniol, cycling vCM, n = 89
    @test MATURATION_SLOPE_E2F6 == 0.396
    @test MATURATION_SLOPE_ECT2 < 0 && MATURATION_SLOPE_E2F6 > 0
    @test MATURATION_SLOPE_ECT2 / MATURATION_SLOPE_E2F6 ≈ -1.4217 atol = 1e-4

    # Inert at M = 0 for any gain -- this is what keeps the axis off by default.
    for g in (0.0, 1.0, 3.0, 10.0)
        @test all(maturation_factors(0.0, g) .== (1.0, 1.0))
    end
    # Ect2 suppressed, E2F6 raised, monotonically in M.
    e_prev, f_prev = 1.0, 1.0
    for M in (0.12, 0.3, 0.5, 0.95)
        e, f = maturation_factors(M, 3.0)
        @test e < e_prev && f > f_prev
        @test e > 0                            # saturating form cannot go negative
        e_prev, f_prev = e, f
    end
    # A linear suppression would have gone negative here; the saturating form does not.
    @test maturation_factors(0.95, 10.0)[1] > 0
end

@testset "maturation closes the abscission arm" begin
    # Tier 1's mechanism: maturation closes E2Fact -> Ect2 -> RhoA -> Midbody. Fate
    # should switch from division to binucleation as M rises, with nothing else changing.
    W4 = (1900.0, 2180.0)
    function at(M)
        en = merge(merge(E2F_SPLIT_ON, CYTOKINESIS_ON), merge(MATURATION_ON, (M = M,)))
        sol, log = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en,
                               record_events = true, cytokinesis = true)
        cyc = classify_cycles(log; window = W4, cytokinesis = true)
        w = findall(x -> W4[1] <= x <= W4[2], sol.t)
        (div = count(c -> c.fate === :Division, cyc),
         bi  = count(c -> c.fate === :Binucleation, cyc),
         ect2 = maximum(sol[species_index("Ect2"), :][w]),
         mb = maximum(sol[species_index("Midbody"), :][w]),
         per = peak_period(sol))
    end
    rs = [at(maturation(c)) for c in context_names()]

    @test issorted([r.ect2 for r in rs], rev = true)   # Ect2 falls with maturation
    @test issorted([r.mb   for r in rs], rev = true)   # and the midbody follows
    @test rs[1].div > 0 && rs[1].bi == 0               # least mature divides
    @test rs[end].div == 0 && rs[end].bi > 0           # adult binucleates
    # The switch happens inside the observed range, not outside it.
    @test any(r -> r.div > 0, rs) && any(r -> r.bi > 0, rs)

    # The arm reads the cycle but does not drive it, so the period must not move with M.
    # Measured spread across the six contexts is 39.3215-39.3422 h, i.e. 0.05 % -- and
    # the two contexts that share M = 0.50 give bit-identical periods, which is what
    # says the residual is step-selection noise rather than a maturation effect.
    pers = [r.per for r in rs]
    @test maximum(pers) - minimum(pers) < 0.1
    @test rs[3].per == rs[4].per        # mouse_p1_invivo and mncm_invitro, both M = 0.50
end

@testset "E2F6 is the exit enforcer, and acts on the period" begin
    # A separate knob from the Ect2 coupling, on a separate observable. Tier 1 labels
    # E2f6 "cell-cycle exit enforcer"; here rising E2F6 lengthens the cycle toward exit.
    function per(ks6)
        en = merge(merge(E2F_SPLIT_ON, CYTOKINESIS_ON),
                   merge(MATURATION_ON, (M = 0.55, ks_E2F6 = ks6)))
        peak_period(solve_tier2(alpha = PUBLISHED_ALPHA, enable = en))
    end
    ps = [per(k) for k in (0.0, 0.005, 0.01, 0.02, 0.05)]
    @test all(isfinite, ps)
    @test issorted(ps)                    # more E2F6, longer cycle
    @test ps[end] > 2 * ps[1]             # and the effect is large
    @test E2F6_EXIT_ON.ks_E2F6 > 0
    @test !haskey(MATURATION_ON, :ks_E2F6)   # the two presets stay separate
end

@testset "M = 0 keeps the reduction bit-exact" begin
    st = 0; worst = 0.0
    for _ in 1:50
        u63, u66 = state(), tier2_state()
        for s in state_names()
            st = (1103515245 * st + 12345) % 2147483648
            v = 0.8 * st / 2147483648
            u63[Symbol(s)] = v; u66[Symbol(s)] = v
        end
        d63, d66 = similar(u63), similar(u66)
        modelDiffEq!(d63, u63, params(), 0.0)
        # Non-zero gain, but M = 0, so the axis must still contribute nothing.
        tier2DiffEq!(d66, u66, tier2_params(maturation_gain = 3.0), 0.0)
        for i in 1:63
            worst = max(worst, abs(d66[i] - d63[i]))
        end
    end
    @test worst == 0.0
end

@testset "REGRESSION: landmark detectors must not chatter" begin
    # A ContinuousCallback on a level a slowly-varying signal lingers near will re-trigger
    # every step. Total geminin spends 12-19 % of the cycle within +/-10 % of the FUCCI
    # cutoff, and with E2F6 repression on that logged 995,356 crossings and drove the
    # solver to MaxIters at t = 1518 (2e6 steps). Fixed with a Schmitt-trigger deadband
    # plus a refractory guard. This pins both.
    thr = EventThresholds()
    @test thr.hysteresis > 0
    @test thr.refractory > 0

    en = merge(merge(E2F_SPLIT_ON, CYTOKINESIS_ON),
               merge(MATURATION_ON, (M = 0.55, ks_E2F6 = 0.05)))
    sol, log = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en,
                           record_events = true, cytokinesis = true)
    @test string(sol.retcode) == "Success"
    @test sol.t[end] == 2500.0
    @test length(sol.t) < 100_000                 # was 1,995,044

    # One of each landmark per cycle, not thousands.
    n = length(log.restriction)
    @test 0 < n < 200
    for v in (log.s_entry, log.neb, log.anaphase, log.mitotic_exit, log.geminin_on)
        @test abs(length(v) - n) <= 2
    end
    # Recorded times still respect the refractory gap.
    for v in (log.restriction, log.geminin_on, log.mitotic_exit)
        length(v) > 1 && @test minimum(diff(sort(v))) > thr.refractory
    end
end

end # testset


# ---------------------------------------------------------------------------
# Phase 2 step 5: oxidative stress -> DDR -> Ccng1 -> the mitotic-entry brake.
# ---------------------------------------------------------------------------

const FULL = merge(merge(merge(E2F_SPLIT_ON, CYTOKINESIS_ON), MATURATION_ON), DDR_ON)
const W5 = (1900.0, 2400.0)

fates_at(over; base = FULL, tspan = (0.0, 2600.0)) = begin
    sol, log = solve_tier2(alpha = PUBLISHED_ALPHA, enable = merge(base, over),
                           tspan = tspan, record_events = true, cytokinesis = true)
    cyc = classify_cycles(log; window = W5, cytokinesis = true)
    (sol = sol, log = log,
     D = count(c -> c.fate === :Division, cyc),
     B = count(c -> c.fate === :Binucleation, cyc),
     P = count(c -> c.fate === :Polyploidization, cyc),
     n_neb = length(trim(log, W5).neb), n_s = length(trim(log, W5).s_entry))
end

@testset "CmTier2 — Phase 2 step 5 (DDR)" begin

@testset "context inputs honour Tier 1's default_on semantics" begin
    # REGRESSION. An input on the manifest's default_on list that a context does not name
    # sits at 1.0, not 0. Reading unlisted-as-zero gave the adult heart no oxidative
    # stress and produced adult cycling FASTER than P0 (39.3 h vs 52.2 h) -- backwards,
    # and exactly the kind of silent mis-mapping that would poison Phase 3.
    @test context_params("adult").ROSenv == 1.0      # not named by `adult`; default_on
    @test context_params("adult").InVitro == 0.0     # named explicitly
    @test context_params("mouse_p1_invivo").ROSenv == 0.20
    @test context_params("hipsc_cm").InVitro == 1.0
    @test context_params("mncm_invitro").M == 0.50
    # InVitro is NOT on default_on, so an unlisted InVitro is 0 rather than 1.
    manifest = contexts()
    @test !haskey(manifest["adult"]["inputs"], "ROSenv")
end

@testset "the arm is inert until switched on, and reduction survives" begin
    @test "Ccng1" in TIER2_SPECIES
    @test species_index("Ccng1") > 63
    sol = solve_tier2(alpha = PUBLISHED_ALPHA)
    @test abs(sol(2000.0)[species_index("Ccng1")]) < 1e-30

    # Bit-exact reduction with the DDR block present and ROS inputs set: the arm is
    # still gated off by ks_ATM_ROS = 0, and fix_p53_massbalance defaults off so the
    # claim holds for ANY state rather than only the reachable ones where Chk2p = 0.
    st = 0; worst = 0.0
    for _ in 1:50
        u63, u66 = state(), tier2_state()
        for s in state_names()
            st = (1103515245 * st + 12345) % 2147483648
            v = 0.8 * st / 2147483648
            u63[Symbol(s)] = v; u66[Symbol(s)] = v
        end
        # Chk2p must be zero for the reduction to hold; the inherited model guarantees
        # that dynamically, so pin it here rather than randomising it.
        d63, d66 = similar(u63), similar(u66)
        modelDiffEq!(d63, u63, params(), 0.0)
        tier2DiffEq!(d66, u66, tier2_params(ROSenv = 0.5, InVitro = 1.0), 0.0)
        for i in 1:63
            worst = max(worst, abs(d66[i] - d63[i]))
        end
    end
    @test worst == 0.0
end

@testset "the p53 mass-balance defect is fixed" begin
    # d.p53p gained p53p at a rate independent of available p53. With the arm live that
    # is no longer harmless. The fix restores the *p53 factor that d.p53's loss term has.
    src = read(joinpath(@__DIR__, "..", "src", "inherited", "diff_eqns.jl"), String)
    @test occursin("d.p53p = ((kf_p53p*Chk2p)-kr_p53p*PPase*p53p) * α", src)  # still broken upstream

    u = tier2_state(); u.Chk2p = 0.5; u.p53 = 0.3; u.PPase = 1.0; u.p53p = 0.1
    d_fixed  = similar(u); d_broken = similar(u)
    tier2DiffEq!(d_fixed,  u, tier2_params(fix_p53_massbalance = 1.0), 0.0)
    tier2DiffEq!(d_broken, u, tier2_params(fix_p53_massbalance = 0.0), 0.0)
    @test d_fixed.p53p != d_broken.p53p

    # Fixed: the p53 -> p53p flux in d.p53p matches the loss term in d.p53 exactly, so
    # the pair is conserved under phosphorylation. Broken: it does not.
    p = tier2_params()
    flux = p.kf_p53p * u.Chk2p * u.p53 * p.α
    @test d_fixed.p53p ≈ (p.kf_p53p * u.Chk2p * u.p53 - p.kr_p53p * u.PPase * u.p53p) * p.α
    @test !isapprox(d_broken.p53p, d_fixed.p53p; rtol = 1e-6)
end

@testset "ROS drives the DDR cascade" begin
    peak(sol, s) = (w = findall(x -> W5[1] <= x <= W5[2], sol.t);
                    maximum(sol[species_index(s), :][w]))
    off = fates_at((M = 0.50, ROSenv = 0.20, InVitro = 1.0, ks_ATM_ROS = 0.0))
    on  = fates_at((M = 0.50, ROSenv = 0.20, InVitro = 1.0))
    for s in ("ATMp", "Chk2p", "p53p", "p21", "Ccng1")
        @test peak(off.sol, s) < 1e-10        # cascade silent without the ROS input
        @test peak(on.sol, s) > 0             # and live with it
    end
    # Monotone in the input, all the way down the cascade.
    lo = fates_at((M = 0.50, ROSenv = 0.20, InVitro = 0.0))
    hi = fates_at((M = 0.50, ROSenv = 0.20, InVitro = 1.0))
    for s in ("ATMp", "p53p", "Ccng1")
        @test peak(hi.sol, s) > peak(lo.sol, s)
    end
end

@testset "culture closes the mitotic-entry brake — the held-out contrast" begin
    # Tier 1's 2x2: maturation closes the abscission arm, culture closes the
    # mitotic-entry brake. Step 4 gave the first half; this is the second.
    #
    # mouse_p1_invivo and mncm_invitro sit at the SAME M = 0.50 and differ only in
    # culture, so this contrast isolates the ROS arm. Both are held out -- Tier 1 fitted
    # neither, and Tier 2 has nothing fitted to either.
    p1   = fates_at(context_params("mouse_p1_invivo"))
    mncm = fates_at(context_params("mncm_invitro"))

    @test p1.B > 0 && p1.P == 0            # in vivo: binucleation
    @test mncm.P > 0 && mncm.B == 0        # in culture: polyploidization
    # The mechanism must be a blocked mitotic entry with S phase intact -- not simply
    # fewer cycles. Polyploidization IS S-without-mitosis.
    @test mncm.n_neb == 0
    @test mncm.n_s > 0
    @test p1.n_neb > 0
end

@testset "the brake acts through MPF, gradedly" begin
    lm(kb) = (r = fates_at((M = 0.50, ROSenv = 0.20, InVitro = 1.0, kf_CCNB_Ccng1 = kb));
              w = findall(x -> W5[1] <= x <= W5[2], r.sol.t);
              (mpf = maximum(r.sol[species_index("CCNB_CDK1"), :][w]),
               lmnap = maximum(r.sol[species_index("LMNAp"), :][w]), neb = r.n_neb))
    rs = [lm(k) for k in (0.0, 10.0, 30.0, 80.0, 200.0)]
    @test issorted([r.mpf for r in rs], rev = true)     # MPF falls with brake strength
    @test issorted([r.lmnap for r in rs], rev = true)   # and NEB follows it down
    @test rs[1].neb > 0 && rs[end].neb == 0
    # NEB stops exactly when LMNAp can no longer reach its threshold.
    @test rs[end].lmnap < EventThresholds().neb
    @test rs[1].lmnap > EventThresholds().neb
end

@testset "KNOWN MISS: hiPSC-CM is predicted polyploid, observed dividing" begin
    # Tier 1 observes Division at hipsc_cm; Tier 2 predicts Polyploidization. Recorded,
    # not tuned away.
    #
    # Diagnosis: hipsc_cm and mncm_invitro carry the SAME oxidative input (ROSenv 0.25,
    # InVitro 1.0), so the DDR brake hits them equally and only M distinguishes them --
    # and M acts on the abscission arm, not on mitotic entry. Tier 1 gets hipsc right
    # because hipsc IS its calibration context: its MitoticEntry gate is fitted to 0.7698
    # there, and MODEL.md says plainly that "the fate layer is fitted exactly and predicts
    # nothing" at that context. Tier 2 has nothing fitted to hipsc, so the miss is a
    # genuine prediction failure and says the model lacks whatever keeps the DDR response
    # weak in immature cardiomyocytes.
    #
    # NOT fixed by gating the DDR arm on !Maturation. That is one free parameter fitted to
    # one outcome, and Tier 1 tested the directly analogous move -- gating the clonidine
    # response on !Maturation -- and REJECTED it: mean fold error went 26 % -> 55 %, worse
    # than doing nothing.
    h = fates_at(context_params("hipsc_cm"))
    @test h.P > 0                    # the miss, asserted so it cannot vanish silently
    @test h.D == 0
    # The two in-vitro contexts really do receive identical oxidative input.
    @test context_params("hipsc_cm").ROSenv == context_params("mncm_invitro").ROSenv
    @test context_params("hipsc_cm").InVitro == context_params("mncm_invitro").InVitro
end

@testset "the Ccng1 brake has a parameter degeneracy worth declaring" begin
    # Only the product ks_Ccng1_p53 * kf_CCNB_Ccng1 / kd_Ccng1 sets the brake strength:
    # Ccng1 is at quasi-steady state, so halving its synthesis and doubling its effect
    # is the same model. Phase 3 must fit ONE effective parameter here, not three.
    a = fates_at((M = 0.50, ROSenv = 0.20, InVitro = 1.0,
                  ks_Ccng1_p53 = 0.80, kf_CCNB_Ccng1 = 200.0))
    b = fates_at((M = 0.50, ROSenv = 0.20, InVitro = 1.0,
                  ks_Ccng1_p53 = 0.40, kf_CCNB_Ccng1 = 400.0))
    @test a.P == b.P && a.B == b.B && a.D == b.D
    @test a.n_neb == b.n_neb
end

end # testset


# ---------------------------------------------------------------------------
# Phase 3: the parameter budget and the population layer.
# ---------------------------------------------------------------------------

@testset "CmTier2 — Phase 3" begin

@testset "the parameter budget is mechanized, not just documented" begin
    # Tier 1 enforces its budget with a linter that fails the build. A budget that is
    # only written down in a README is a budget that drifts.
    @test isempty(lint_budget())

    # Completeness: every Tier-2 parameter has a declared provenance, and nothing is
    # declared that does not exist. This is what stops the budget drifting silently when
    # a later step adds a parameter.
    declared = Set(keys(PARAM_PROVENANCE))
    actual = Set(setdiff(collect(keys(tier2_params())), collect(keys(params()))))
    @test declared == actual

    # The count the paper has to defend.
    @test total_fitted() == 6
    @test total_fitted() <= FITTED_BUDGET
    @test FITTED_BUDGET < FITTED_CEILING      # headroom is the argument
    @test FITTED_CEILING == 13                # MODEL.md 3.5

    # The linter must actually fail when the budget is exceeded, or it proves nothing.
    @test !isempty(lint_budget(budget = 2))
    # Untying the E2F repressors costs a parameter, and the linter must see that.
    @test length(lint_budget(budget = total_fitted(), tied = false)) > 0

    # The FUCCI cutoff is fitted but is not a tier2_params entry, so it can only be
    # caught by being declared separately. Assert it is.
    @test :FUCCI_THRESHOLD in FITTED_OBSERVATION_PARAMS

    # Every fitted parameter carries a justification, not just a label.
    for k in fitted_params()
        @test length(PARAM_PROVENANCE[k][2]) > 20
    end
    @test occursin("STRUCTURAL", budget_report())    # the grey zone is stated
end

@testset "heterogeneity perturbs abundances, not chemistry" begin
    p0 = tier2_params()
    p = heterogeneous_params(p0, 0.2, Random.Xoshiro(1))
    nchanged = count(k -> p[k] != p0[k], keys(p0))
    @test nchanged > 10
    # Only synthesis rates move. A rate constant is a property of the molecule; an
    # abundance is a property of the cell.
    for k in keys(p0)
        startswith(String(k), "ks_") || @test p[k] == p0[k]
    end
    # Lognormal: strictly positive, median 1, so no draw can make synthesis negative.
    draws = [heterogeneous_params(p0, 0.5, Random.Xoshiro(i))[:ks_CCNE_E2F] for i in 1:400]
    @test all(>(0), draws)
    @test 0.7 < Statistics.median(draws) / p0[:ks_CCNE_E2F] < 1.4
    # sigma = 0 is exactly the unperturbed model.
    @test heterogeneous_params(p0, 0.0, Random.Xoshiro(1)) == p0
end

@testset "ensemble sizing is derived from the experiment's own error" begin
    # There is no value in driving Monte-Carlo error far below the uncertainty in the
    # measurement being compared against. Murganti's rarest fate is 1.40 % of 570 cells.
    @test required_n(0.0140) == 57_000       # TODO.md item 4's "1e4-1e5"
    @test required_n(0.0140, rel = 1.0) == 570   # matching the experiment exactly

    # It does NOT depend on p: both errors are binomial in the same p, so it cancels and
    # n = n_obs/rel^2. One ensemble size buys the same error ratio for every fate at
    # once, rather than having to be sized against the rarest.
    @test required_n(0.05) == required_n(0.0140) == required_n(0.9035)
    @test required_n(0.0140, rel = 0.2) == required_n(0.0140) / 4
end

@testset "the ensemble runs, and its noise scale is measured not chosen" begin
    full = merge(merge(merge(E2F_SPLIT_ON, CYTOKINESIS_ON), MATURATION_ON), DDR_ON)
    en = merge(full, context_params("mouse_p1_invivo"))

    # Deterministic limit: sigma = 0 gives one fate and zero spread.
    r0 = run_ensemble(n = 24, sigma = 0.0, enable = en)
    @test r0.failed == 0
    @test r0.duration_cv ≈ 0.0 atol = 1e-9
    @test maximum(values(r0.fractions)) == 1.0        # every cell agrees
    @test sum(values(r0.fractions)) ≈ 1.0

    # Heterogeneity produces GRADED fractions -- the capability Tier 1 structurally
    # cannot have, since its fates are a product of steady-state gate activities.
    r1 = run_ensemble(n = 96, sigma = 0.10, enable = en)
    @test sum(values(r1.fractions)) ≈ 1.0
    @test count(v -> v > 0, values(r1.fractions)) >= 2
    @test r1.duration_cv > r0.duration_cv

    # The width is pinned by Baniol's measured duration CV, NOT by the fate fractions --
    # those are the hold-out, and fitting the noise to them would make the headline
    # prediction circular. Measured: CV rises 0.000 -> 0.155 -> 0.316 over
    # sigma 0.00 -> 0.01 -> 0.02, so the target 0.265 lands near 0.016.
    @test MEASURED_DURATION_CV ≈ 4.0 / 15.1
    @test 0.25 < MEASURED_DURATION_CV < 0.28
    cvs = [run_ensemble(n = 96, sigma = s, enable = en).duration_cv
           for s in (0.0, 0.01, 0.02)]
    @test issorted(cvs)
    @test cvs[3] > MEASURED_DURATION_CV > cvs[2]   # the target is bracketed

    # Reproducible: same seed, same answer, independent of thread scheduling.
    a = run_ensemble(n = 48, sigma = 0.05, enable = en, seed = 7)
    b = run_ensemble(n = 48, sigma = 0.05, enable = en, seed = 7)
    @test a.fractions == b.fractions
    @test run_ensemble(n = 48, sigma = 0.05, enable = en, seed = 8).counts != a.counts
end

@testset "pathological draws fail fast and are reported, not hidden" begin
    # Heterogeneity draws occasional parameter sets that are very stiff, and
    # Threads.@threads makes wall-clock the SLOWEST cell -- at a 2e6 step cap one bad
    # draw stalled a 64-cell run past 13 minutes. The cap is sized from measurement: a
    # healthy cell needs ~7,600 steps.
    @test MAXITERS_PER_CELL == 100_000
    @test MAXITERS_PER_CELL > 10 * 7_600
    full = merge(merge(merge(E2F_SPLIT_ON, CYTOKINESIS_ON), MATURATION_ON), DDR_ON)
    r = run_ensemble(n = 96, sigma = 0.10,
                     enable = merge(full, context_params("mouse_p1_invivo")))
    @test r.scored + r.failed == r.n          # nothing silently dropped
    @test r.failed / r.n < 0.05               # and failures stay rare
    # Fractions are over SCORED cells, so a failure cannot masquerade as a fate.
    @test sum(values(r.fractions)) ≈ 1.0
end

end # testset
