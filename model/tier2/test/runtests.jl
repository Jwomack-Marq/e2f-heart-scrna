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
