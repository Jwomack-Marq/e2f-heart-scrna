# Named biological contexts, read from Tier 1's manifest.
#
# Deliberately NOT duplicated here. `../cmcycle/data/cmfate_model.toml` is the single
# definition of what "P7 mouse" means, and both tiers read it. If a context drifts
# between tiers that is a bug, not a configuration — and cross-tier comparison is only
# meaningful if the two are answering about the same cell.

"""Path to Tier 1's manifest, the shared source of context definitions."""
const CMFATE_MANIFEST = normpath(joinpath(@__DIR__, "..", "..",
                                          "cmcycle", "data", "cmfate_model.toml"))

"""
    contexts() -> Dict{String,Any}

The six named contexts from Tier 1's manifest, each a label plus an input vector.

Tier 2 currently consumes only `Maturation`; the remaining inputs (`InVitro`,
`MechLoad`, `BetaAR`, `ROSenv`, `Nrg1`, `IGF1`) belong to the signalling layer that
arrives with step 5, and are carried through untouched so nothing is silently dropped.

| context | M |
|---|---|
| `hipsc_cm` | 0.12 |
| `mouse_p0_invivo` | 0.30 |
| `mncm_invitro` | 0.50 |
| `mouse_p1_invivo` | 0.50 |
| `mouse_p7_invivo` | 0.55 |
| `adult` | 0.95 |
"""
function contexts()
    isfile(CMFATE_MANIFEST) ||
        error("Tier 1 manifest not found at $(CMFATE_MANIFEST). Tier 2 reads context " *
              "definitions from Tier 1 rather than duplicating them.")
    return TOML.parsefile(CMFATE_MANIFEST)["contexts"]
end

"""
    maturation(context) -> Float64

The maturation coordinate `M` for a named context, from Tier 1's manifest.

Contexts that do not name `Maturation` inherit Tier 1's `default_on` behaviour, where an
unlisted default-on input sits at 1.0 — but all six name it explicitly, so this is a
guard rather than a code path in use.
"""
function maturation(context::AbstractString)
    ctx = contexts()
    haskey(ctx, context) ||
        error("unknown context $(repr(context)); known: $(sort(collect(keys(ctx))))")
    return Float64(get(ctx[context]["inputs"], "Maturation", 1.0))
end

"""
    context_names() -> Vector{String}

Context names sorted by maturation, which is the axis Tier 2 varies.
"""
context_names() = sort(collect(keys(contexts())); by = maturation)
