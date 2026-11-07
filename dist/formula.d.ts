/**
 * Composition → formula: turn the elemental analysis of a compound (mass
 * percentages, or plain masses) into an empirical formula, and — when a molar
 * mass is known — the molecular formula.
 *
 * The arithmetic is the classic divide-by-atomic-mass / divide-by-smallest /
 * scale-to-integers loop, which agents get wrong in predictable ways: they
 * forget that the smallest ratio itself may not be 1, they round 1.5 to 2
 * without noticing, and they happily report a "molecular formula" from a molar
 * mass that is not an integer multiple of the empirical mass. This module
 * reports the multiplier, the residual and an explicit `exact` flag instead of
 * quietly rounding.
 *
 * Deterministic rules (documented in the README):
 *   - `percent` is mass percent and must sum to 100 ± 1 (a typo is an error,
 *     not something to normalise away silently);
 *   - subscripts come from the SMALLEST multiplier k ∈ [1, 8] whose scaled
 *     ratios are all within ±0.05 of an integer;
 *   - if no k fits, the closest fit is reported with `exact: false`;
 *   - formulas are written in Hill order (C, H, then alphabetical; alphabetical
 *     when there is no carbon) so the same analysis always prints the same
 *     string regardless of input order.
 *
 * @module dsh-chem/formula
 */
/** Largest number of elements accepted in one composition. */
export declare const MAX_ELEMENTS = 12;
/** Largest integer multiplier searched for (empirical subscripts are coprime). */
export declare const MAX_MULTIPLIER = 8;
/** |ratio × k − round(ratio × k)| must be within this to count as an integer. */
export declare const RATIO_TOLERANCE = 0.05;
/** A molar mass within this relative distance of an integer multiple is accepted. */
export declare const MOLECULAR_TOLERANCE = 0.02;
export interface CompositionPart {
    /** Element symbol, e.g. "C" (single-letter lowercase is corrected, as elsewhere). */
    symbol: string;
    /** Mass percent (mode "percent") or mass in grams (mode "masses"). */
    value: number;
}
export interface CompositionInput {
    /** "percent" = the values are mass percentages (must sum to 100 ± 1); "masses" = grams. */
    mode: 'percent' | 'masses';
    parts: CompositionPart[];
    /** Measured molar mass in g/mol — turns the answer into a molecular formula. */
    molarMass?: number;
}
export interface CompositionElementOut {
    symbol: string;
    /** Mass used in the calculation (grams, or grams per 100 g for `percent`). */
    massG: number;
    atomicMass: number;
    moles: number;
    /** moles ÷ smallest moles (before scaling). */
    ratio: number;
    /** Integer subscript in the empirical formula. */
    subscript: number;
    /** Mass percent of this element in the empirical formula (round-trip check). */
    massFraction: number;
}
export type CompositionResult = {
    ok: true;
    mode: 'percent' | 'masses';
    empiricalFormula: string;
    empiricalMolarMass: number;
    multiplier: number;
    maxDeviation: number;
    exact: boolean;
    elements: CompositionElementOut[];
    molecularFormula?: string;
    molecularMultiplier?: number;
    molecularExact?: boolean;
    measuredMolarMass?: number;
    impliedMolarMass?: number;
    molarMassErrorPercent?: number;
    note?: string;
} | {
    ok: false;
    error: string;
};
/** Hill-system ordering: C, then H, then alphabetical; alphabetical without C. */
export declare function hillOrder(symbols: string[]): string[];
/** Render subscripts as a formula string in Hill order (subscript 1 omitted). */
export declare function renderFormula(counts: Record<string, number>): string;
/**
 * Resolve an elemental analysis → empirical (and molecular) formula.
 *
 * `mode` says what the numbers mean ("percent" = mass %, "masses" = grams);
 * `molarMass` is optional.
 */
export declare function compositionToFormula(input: CompositionInput): CompositionResult;
//# sourceMappingURL=formula.d.ts.map