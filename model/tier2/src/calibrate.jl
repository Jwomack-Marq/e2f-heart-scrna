# Joint calibration of the fitted parameters.
#
# ## Why joint, and not one at a time
#
# The five fitted parameters were originally placed by hand, one after another, and that
# method demonstrably does not converge on this system. Twice an upstream change
# invalidated a downstream estimate:
#
#   * making CPC recruitment switch-like raised the whole cytokinesis arm, which pushed
#     the maturation switch outside the observed range and stale'd `maturation_gain`;
#   * re-estimating `maturation_gain` then moved the `ks_Ect2_E2F` window from
#     [0.70, 0.85] to [0.85, 1.15].
#
# The parameters are coupled through the same three observables, so they have to be
# estimated together against one loss.
#
# ## What is fitted against, and what is held out
#
# Calibration set — deliberately DISJOINT from Tier 1's four fitted fate fractions:
#
#   * corrected in-vivo cycling fractions, P0 17.7 % and P7 5.5 % (`../TODO.md` item 6;
#     the correction applies Baniol's own Suppl 1G Ki-67 co-staining, and recovering the
#     17.7 % figure is what gave Tier 1 its strongest unfitted validation)
#   * FUCCI-timed S/G2/M duration, Baniol P0 mouse 15.1 h
#   * the regeneration window: P0 cardiomyocytes divide, P7 largely do not
#
# HELD OUT: Murganti Fig 2A's four fate fractions. Nothing here may touch them, or the
# headline prediction becomes a restatement.

"""Corrected in-vivo cycling fractions (TODO.md item 6, via Baniol Suppl 1G)."""
const TARGET_CYCLING = (p0 = 0.177, p7 = 0.055)

"""FUCCI-timed S/G2/M, Baniol P0 mouse (15.1 +/- 4.0 h)."""
const TARGET_SG2M_P0 = 15.1

"""
    calibration_loss(theta; n, sigma, seed) -> NamedTuple

Scalar loss for one parameter vector, plus its components so a fit can be read rather
than merely trusted.

`theta` is a NamedTuple over the fitted parameters. Terms are scaled by their own target
so that a fraction and a duration contribute comparably; the regeneration-window term is
a one-sided hinge, because "P0 divides more than P7" is an inequality rather than a
number to hit.
"""
function calibration_loss(theta::NamedTuple; n::Int = 200, sigma::Real = 0.016,
                          seed::Int = 20260816)
    base = merge(CALIBRATED, theta)
    # ks_E2F8_E2F is tied to ks_E2F7_E2F: Baniol measure the two at +0.46 and +0.49, so
    # resolving them separately would spend a parameter the data cannot support.
    haskey(theta, :ks_E2F7_E2F) && (base = merge(base, (ks_E2F8_E2F = theta.ks_E2F7_E2F,)))

    out = Dict{Symbol,Float64}()
    for (key, ctx) in ((:p0, "mouse_p0_invivo"), (:p7, "mouse_p7_invivo"))
        r = run_ensemble(n = n, sigma = sigma, seed = seed,
                         enable = merge(base, context_params(ctx)))
        out[Symbol(key, :_cyc)] = 1 - r.fractions[:Quiescent]
        out[Symbol(key, :_div)] = r.fractions[:Division]
        out[Symbol(key, :_dur)] = r.duration_mean
        out[Symbol(key, :_fail)] = r.failed / max(r.n, 1)
    end

    l_cyc = ((out[:p0_cyc] - TARGET_CYCLING.p0) / TARGET_CYCLING.p0)^2 +
            ((out[:p7_cyc] - TARGET_CYCLING.p7) / TARGET_CYCLING.p7)^2
    l_dur = isfinite(out[:p0_dur]) ?
            ((out[:p0_dur] - TARGET_SG2M_P0) / TARGET_SG2M_P0)^2 : 4.0
    # Regeneration window: P0 should divide more than P7. One-sided.
    l_win = max(0.0, 0.10 - (out[:p0_div] - out[:p7_div]))^2 / 0.10^2
    # Failed integrations are a modelling problem, not a free lunch.
    l_fail = 4 * (out[:p0_fail] + out[:p7_fail])

    total = l_cyc + l_dur + l_win + l_fail
    return (total = total, cyc = l_cyc, dur = l_dur, win = l_win, fail = l_fail,
            detail = out, theta = theta)
end

"""
Search grids for [`calibrate_joint`](@ref).

Ranges bracket the hand-placed values rather than starting from them, so the fit is free
to reject them. `ks_E2F6` runs down to 0 because whether cell-cycle exit is needed at all
is itself a question the cycling fractions should answer.
"""
const DEFAULT_GRIDS = (
    ks_E2F6         = (0.0, 0.005, 0.01, 0.02, 0.04, 0.08),
    maturation_gain = (3.0, 5.0, 8.0, 12.0),
    ks_Ect2_E2F     = (0.6, 0.85, 1.0, 1.3, 1.8),
    kf_CCNB_Ccng1   = (50.0, 120.0, 200.0, 350.0),
    ks_E2F7_E2F     = (0.05, 0.10, 0.20, 0.40),
)

"""
    calibrate_joint(; grids, passes, n, sigma, verbose) -> NamedTuple

Coordinate descent over the fitted parameters.

Coordinate descent rather than a gradient method because each evaluation is an ensemble
of stochastic ODE solves — there is no usable gradient — and rather than a full grid
because five parameters at five levels is 3,125 evaluations at ~20 s each.

It is a local method and makes no claim to a global optimum. What it does provide, and
what hand-placement did not, is that every parameter is re-examined after every other one
moves, so the cascade that stale'd two estimates cannot recur silently.
"""
function calibrate_joint(; grids = DEFAULT_GRIDS, passes::Int = 2, n::Int = 200,
                           sigma::Real = 0.016, verbose::Bool = true)
    theta = NamedTuple()
    best = calibration_loss(theta; n = n, sigma = sigma)
    verbose && @info "start" loss=round(best.total, digits=4)
    history = [(pass = 0, param = :start, value = NaN, loss = best.total)]

    for pass in 1:passes, (param, values) in pairs(grids)
        cand = best
        for v in values
            trial = calibration_loss(merge(best.theta, (; param => v)); n = n, sigma = sigma)
            trial.total < cand.total && (cand = trial)
        end
        if cand.total < best.total
            best = cand
            push!(history, (pass = pass, param = param,
                            value = get(best.theta, param, NaN), loss = best.total))
            verbose && @info "improved" pass param value=get(best.theta, param, NaN) loss=round(best.total, digits=4)
        end
    end
    return (theta = best.theta, loss = best, history = history)
end

