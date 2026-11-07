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
export type EquilibriumMode = 'weak_acid' | 'weak_base';
export interface EquilibriumInput {
    mode: EquilibriumMode;
    /** Ka for a weak acid, Kb for a weak base (mol/L). */
    k: number;
    /** Formal concentration of the acid/base, mol/L. */
    concentration: number;
}
export interface EquilibriumApproximation {
    /** x from the √(K·C) shortcut. */
    x: number;
    /** How much the shortcut overstates the exact x, in percent. */
    errorPercent: number;
}
export type EquilibriumResult = {
    ok: true;
    mode: EquilibriumMode;
    constantK: number;
    pK: number;
    concentration: number;
    /** Amount dissociated, mol/L (the "x" of the ICE table). */
    x: number;
    hConc: number;
    ohConc: number;
    pH: number;
    pOH: number;
    degreeOfDissociation: number;
    percentDissociation: number;
    approximation: EquilibriumApproximation;
    /** K recomputed from the answer: x² ÷ (C − x). Must equal `constantK`. */
    kFromResult: number;
    note?: string;
} | {
    ok: false;
    error: string;
};
/** Threshold above which the √(K·C) shortcut is flagged as misleading. */
export declare const APPROXIMATION_WARN_PERCENT = 5;
/**
 * Solve `x² + K·x − K·C = 0` for its positive root, in the stable form.
 */
export declare function dissociatedAmount(k: number, concentration: number): number;
/** Full ICE answer for a weak acid or weak base. */
export declare function weakEquilibrium(input: EquilibriumInput): EquilibriumResult;
//# sourceMappingURL=equil.d.ts.map