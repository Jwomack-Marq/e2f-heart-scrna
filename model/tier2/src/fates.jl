# Fate classification and nuclei/ploidy bookkeeping.
#
# The architectural difference from Tier 1. There, the four fates are a complementary
# product over three steady-state gate ACTIVITIES — a population-average partition that
# sums to 1 by construction. Here a fate is a property of one cell's TRAJECTORY: which
# landmarks it actually reached, in what order. The partition is exact because the
# labels are mutually exclusive events, not because of an algebraic identity.
#
# ## Phase 1 scope, stated plainly
#
# Three of Tier 1's four fates are decidable from the inherited model:
#
#   Quiescent          no S entry
#   Polyploidization   S entry, no nuclear-envelope breakdown  (re-replication)
#   MitoticCompletion  S entry -> NEB -> anaphase -> mitotic exit
#
# `MitoticCompletion` is Tier 1's Division + Binucleation. Splitting it needs the
# abscission decision, which needs Ect2/RhoA/Centralspindlin/AurKB/Anillin/Midbody —
# absent from the inherited model, added in Phase 2. It is deliberately NOT guessed at
# here: a cell that completes mitosis without abscission is binucleate, and nothing in
# the current equations represents the furrow, so there is no honest way to call it.

"""
    CellFate

`:Quiescent`, `:Polyploidization`, or `:MitoticCompletion`. Phase 2 splits the last into
`:Division` and `:Binucleation` on the midbody-resolution criterion.
"""
const FATES_PHASE1 = (:Quiescent, :Polyploidization, :MitoticCompletion)

"""
    Cycle

One traversal, from S-phase entry to its outcome.

Two S/G2/M durations are carried, and the distinction matters for calibration:

  * `sg2m_cyclin` — CycE-CDK2 crossing to mitotic exit. A molecular definition.
  * `sg2m_fucci`  — total geminin crossing the FUCCI cutoff to mitotic exit. **This is
    the comparator for Murganti Fig 2E** (16.38 h dividing, 17.29 h binucleating,
    24.50 h polyploidising), because their measurement is a FUCCI trace: mAG (geminin)
    appearing, to the division event. Scoring a cyclin-defined duration against a
    reporter-defined measurement would be a category error, and it is not a small one —
    the two differ by ~8 h in the inherited model.

Tier 1 can produce neither, its `tau` being a relaxation constant rather than a duration.
"""
struct Cycle
    fate::Symbol
    s_entry::Float64
    geminin_on::Union{Float64,Nothing}
    neb::Union{Float64,Nothing}
    anaphase::Union{Float64,Nothing}
    exit::Union{Float64,Nothing}
    sg2m_cyclin::Union{Float64,Nothing}  # CycE crossing -> mitotic exit
    sg2m_fucci::Union{Float64,Nothing}   # geminin appearing -> mitotic exit
    mitosis::Union{Float64,Nothing}      # NEB -> mitotic exit
end

"""
    classify_cycles(log; window) -> Vector{Cycle}

Reconstruct cycles from an [`EventLog`](@ref) and label each one.

## Cycles open at the restriction point, not at CycE activity

The obvious choice — open a cycle when CycE-CDK2 crosses — is **leaky**, and measurably
so. Under strong CDK4/6 inhibition (`con_ABE = 10`) `ppRB` peaks at 0.013, far below the
restriction threshold, so the cell never commits; yet CycE-CDK2 still reaches 0.52 and
crosses a 0.299 marker, because `d.CCNE` carries a `ks_CCNE_pRBE2F` term and a basal
synthesis rate that do not require Rb hyperphosphorylation. A CycE-only rule therefore
reports S-phase entry in a cell that is arrested in G1, which is the one error a fate
classifier must not make.

So a cycle opens at an `ppRB` up-crossing — passage of the restriction point, the
commitment step — and `s_entry` is the first CycE crossing after it. Cells that
accumulate CycE without committing produce no cycles and read as quiescent, which is
what they are.

Landmarks are attributed to the traversal they fall inside rather than to whichever came
nearest in time. A cycle that opens but does not close before the window ends is
dropped: it is right-censored, and counting it would bias durations downward exactly the
way Tier 1's `hipsc_cm` reservoir check describes (the reported 24.5 h there is a
right-censored lower bound, not a mean).
"""
function classify_cycles(log::EventLog; window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing)
    l = window === nothing ? log : trim(log, window)
    commits = sort(l.restriction)
    length(commits) < 2 && return Cycle[]

    first_in(v, lo, hi) = (i = findfirst(t -> lo < t <= hi, v); i === nothing ? nothing : v[i])

    cycles = Cycle[]
    for k in 1:(length(commits) - 1)
        t0, t1 = commits[k], commits[k+1]
        sen  = first_in(sort(l.s_entry),      t0, t1)
        gem  = first_in(sort(l.geminin_on),   t0, t1)
        neb  = first_in(sort(l.neb),          t0, t1)
        ana  = first_in(sort(l.anaphase),     t0, t1)
        ex   = first_in(sort(l.mitotic_exit), t0, t1)

        # No S entry after committing: the cycle stalled before replication. Not a fate
        # in Tier 1's partition, so it is not forced into one.
        sen === nothing && continue

        fate = neb === nothing ? :Polyploidization : :MitoticCompletion
        sg2m_c = ex === nothing ? nothing : ex - sen
        sg2m_f = (gem === nothing || ex === nothing) ? nothing : ex - gem
        mit    = (neb === nothing || ex === nothing) ? nothing : ex - neb
        push!(cycles, Cycle(fate, sen, gem, neb, ana, ex, sg2m_c, sg2m_f, mit))
    end
    return cycles
end

"""
    quiescent(log; window) -> Bool

True when the restriction point was never passed in `window` — the cell never committed,
so it has no cycles to classify.

This is the direct counterpart of Tier 1's `!SPhase => Quiescent` rule, and it is keyed
on `restriction` rather than `s_entry` for the reason given in [`classify_cycles`](@ref):
CycE-CDK2 accumulates under CDK4/6 inhibition without commitment, so an `s_entry`-keyed
test calls an arrested cell cycling.

Kept separate from [`classify_cycles`](@ref) because quiescence is the *absence* of the
event that opens a cycle, so it cannot be a row in a per-cycle table.
"""
quiescent(log::EventLog; window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing) =
    isempty((window === nothing ? log : trim(log, window)).restriction)

"""
    Bookkeeping

Nuclei, DNA content in C, and cell count — `../TODO.md` item 3's
"nuclei-and-ploidy bookkeeping".

Tier 1 has no state that can carry these, which is why it cannot be scored against the
two `ploidy` rows (mNCM tetraploid 26.52 % -> 38.93 % on clonidine) or the two
`nucleation` rows (P7 mononucleated 22.3 % -> 14.8 %) in
`../cmcycle/data/cmcycle_targets.csv`.
"""
struct Bookkeeping
    cells::Float64
    nuclei::Float64      # per original cell
    dna_content::Float64 # in C; 2C is unreplicated diploid
end
Bookkeeping() = Bookkeeping(1.0, 1.0, 2.0)

"""
    bookkeep(cycles; assume_abscission=true) -> Bookkeeping

Walk a fate sequence and accumulate cell, nucleus and DNA counts.

S-phase entry doubles DNA content. A `:MitoticCompletion` then either divides the cell
(abscission succeeds) or leaves one cell with two nuclei (it fails); `:Polyploidization`
leaves the replicated DNA in a single nucleus, so C keeps doubling.

`assume_abscission` is the Phase 1 placeholder for the decision Phase 2 makes
mechanically from the midbody. It is a keyword rather than a hidden default so that
every call site has to state which convention it is using, and so the Phase 2 diff is
visible.
"""
function bookkeep(cycles::AbstractVector{Cycle}; assume_abscission::Bool = true)
    b = Bookkeeping()
    cells, nuclei, dna = b.cells, b.nuclei, b.dna_content
    for c in cycles
        dna *= 2                                    # S phase
        if c.fate === :MitoticCompletion
            if assume_abscission
                cells *= 2
                dna /= 2                            # each daughter gets 2C
            else
                nuclei += 1                         # binucleate, DNA split between nuclei
                dna /= 2
            end
        end
        # :Polyploidization keeps the doubled DNA in one nucleus
    end
    return Bookkeeping(cells, nuclei, dna)
end

"""
    fate_summary(log; window, assume_abscission=true) -> Dict

Fate counts, fractions, mean durations and final bookkeeping for one trajectory.

The fractions are over completed cycles. For a single deterministic trajectory they are
0 or 1 for whichever fate the limit cycle settles into — the graded fractions Murganti
Fig 2A reports are a *population* quantity and need the Phase 3 ensemble, since
`../TODO.md` item 4 notes fate fractions are tail statistics that a mean-field treatment
cannot reach.
"""
function fate_summary(log::EventLog;
                      window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
                      assume_abscission::Bool = true)
    if quiescent(log; window = window)
        return Dict(:quiescent => true, :n_cycles => 0,
                    :counts => Dict(f => 0 for f in FATES_PHASE1),
                    :book => Bookkeeping())
    end
    cycles = classify_cycles(log; window = window)
    counts = Dict(f => count(c -> c.fate === f, cycles) for f in FATES_PHASE1)
    dur(sel) = (v = [x for x in sel if x !== nothing]; isempty(v) ? nothing : sum(v)/length(v))
    return Dict(
        :quiescent => false,
        :n_cycles  => length(cycles),
        :counts    => counts,
        :fractions => Dict(f => (isempty(cycles) ? 0.0 : counts[f]/length(cycles))
                           for f in FATES_PHASE1),
        :mean_sg2m_cyclin => dur(c.sg2m_cyclin for c in cycles),
        :mean_sg2m_fucci  => dur(c.sg2m_fucci  for c in cycles),
        :mean_mitosis     => dur(c.mitosis     for c in cycles),
        :book => bookkeep(cycles; assume_abscission = assume_abscission),
    )
end
