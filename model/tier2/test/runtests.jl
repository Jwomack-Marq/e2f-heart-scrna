using Test
using CmTier2
using Statistics

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
                  nothing, nothing, nothing) for _ in 1:3]
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
    @test length(tier2_params()) == 218 + 13
    @test tier2_state_names()[1:63] == state_names()
    # New components MUST come last: diff_eqns.jl destructures positionally.
    @test tier2_state_names()[64:end] == collect(TIER2_SPECIES)
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
    function fit_err(en)
        sol = solve_tier2(alpha = PUBLISHED_ALPHA, enable = en)
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
    @test off < 0.01                                        # published is near-perfect
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
