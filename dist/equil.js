/**
 * Weak acid / weak base equilibrium (ICE arithmetic) at 25 °C.
 *
 * The one calculation everyone skips: a weak acid does not dissociate
 * completely, so [H+] is the positive root of x² + Kₐ·x − Kₐ·C = 0 — not
 * √(Kₐ·C). The square-root shortcut is only an approximation that *overstates*
 * x, and the error grows as the solution gets more dilute (6.9 % at 0.001 M
 * ammonia with K_b = 1.8e-5, 0.21 % at 1 M acetic acid). This module solves the
 * quadratic exactly and reports both numbers, so the caller can see whether the
 * shortcut was good enough instead of assuming it was.
 *
 * Numerical note: the root is evaluated in the cancellation-free form
 * x = 2·K·C / (K + √(K² + 4·K·C)) instead of (−K + √(K² + 4·K·C))/2 — the two
 * are algebraically identical, but the latter loses precision when K ≪ C
 * (which is exactly the weak-acid regime).
 *
 * Conventions: Ka given as the acid dissociation constant (Kb for bases), C as
 * the formal (analytical) concentration of the acid/base in mol/L, water
 * autoionization Kw = 1e-14 at 25 °C. The tool does not model activity
 * coefficients, ionic strength, polyprotic acids or temperature dependence.
 *
 * @module dsh-chem/equil
 */
import { KW_25C, roundOut } from "./core.js";
/** Threshold above which the √(K·C) shortcut is flagged as misleading. */
export const APPROXIMATION_WARN_PERCENT = 5;
/**
 * Solve `x² + K·x − K·C = 0` for its positive root, in the stable form.
 */
export function dissociatedAmount(k, concentration) {
    const disc = k * k + 4 * k * concentration;
    return (2 * k * concentration) / (k + Math.sqrt(disc));
}
/** Full ICE answer for a weak acid or weak base. */
export function weakEquilibrium(input) {
    const { mode, k, concentration } = input;
    if (mode !== 'weak_acid' && mode !== 'weak_base') {
        return { ok: false, error: `unknown mode "${mode}" (expected weak_acid or weak_base)` };
    }
    if (!Number.isFinite(k))
        return { ok: false, error: `k must be a finite number (got ${k})` };
    if (k <= 0)
        return { ok: false, error: `k must be positive (got ${k})` };
    if (!Number.isFinite(concentration)) {
        return { ok: false, error: `concentration must be a finite number (got ${concentration})` };
    }
    if (concentration <= 0)
        return { ok: false, error: `concentration must be positive (got ${concentration})` };
    const x = dissociatedAmount(k, concentration);
    const isAcid = mode === 'weak_acid';
    const hConc = isAcid ? x : KW_25C / x;
    const ohConc = isAcid ? KW_25C / x : x;
    const pH = -Math.log10(hConc);
    const pOH = -Math.log10(ohConc);
    const approximate = Math.sqrt(k * concentration);
    const approximation = {
        x: roundOut(approximate),
        errorPercent: roundOut(((approximate - x) / x) * 100),
    };
    const result = {
        ok: true,
        mode,
        constantK: roundOut(k),
        pK: roundOut(-Math.log10(k)),
        concentration: roundOut(concentration),
        x: roundOut(x),
        hConc: roundOut(hConc),
        ohConc: roundOut(ohConc),
        pH: roundOut(pH),
        pOH: roundOut(pOH),
        degreeOfDissociation: roundOut(x / concentration),
        percentDissociation: roundOut((x / concentration) * 100),
        approximation,
        kFromResult: roundOut((x * x) / (concentration - x)),
    };
    if (Math.abs(approximation.errorPercent) >= APPROXIMATION_WARN_PERCENT) {
        result.note = `the √(K·C) shortcut gives ${roundOut(approximate, 6)} instead of ${roundOut(x, 6)} — `
            + `${roundOut(approximation.errorPercent, 4)} % too high at this concentration; use the quadratic answer (x)`;
    }
    return result;
}
//# sourceMappingURL=equil.js.map