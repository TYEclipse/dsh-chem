/**
 * Pure chemistry math for dsh-chem: formula parsing, molar mass, mole/particle
 * conversions, solution concentration and pH — zero dependencies, no I/O.
 *
 * All rules are deterministic and documented; agents should call these tools
 * instead of doing molar-mass / dilution / pH arithmetic in their head.
 *
 * @module dsh-chem/core
 */
/** SI Avogadro constant (exact, mol^-1, CODATA 2018). */
export declare const AVOGADRO = 6.02214076e+23;
/** Water autoionization constant at 25 °C (mol²/L²). */
export declare const KW_25C = 1e-14;
export interface ElementCount {
    symbol: string;
    count: number;
}
export interface ParsedGroup {
    /** Flat per-element counts of this group (symbols in first-appearance order). */
    elements: ElementCount[];
    /** Canonical reconstruction of the group, e.g. "(NH4)2" or "Fe". */
    source: string;
}
export type FormulaParse = {
    ok: true;
    /** Main formula groups in source order (non-parenthesized atoms included). */
    groups: ParsedGroup[];
    /** Flat total per-element counts across all groups. */
    totals: ElementCount[];
    /** Hydrate water-group counts if a hydrate separator was present, else null. */
    hydrate: ElementCount[] | null;
    /** Leading integer coefficient if present, else null. */
    coefficient: number | null;
    /** Ionic charge suffix ("2-", "3+", "+") if present, else null. */
    charge: string | null;
    /** Canonical normalized formula (coefficient and charge preserved). */
    normalized: string;
} | {
    ok: false;
    error: string;
};
/**
 * Parse a chemical formula: nested parentheses, hydrate dot (CuSO4·5H2O),
 * leading coefficient (2H2O) and trailing ionic charge (Fe3+, SO4^2-, OH-).
 */
export declare function parseFormula(input: string): FormulaParse;
export interface MassBreakdownEntry {
    symbol: string;
    name: string;
    count: number;
    atomicMass: number;
    /** Absolute mass contribution in g/mol. */
    massContribution: number;
    /** Mass fraction of the total (0–1). */
    massFraction: number;
}
export interface MassResult {
    /** Total molar mass in g/mol. */
    molarMass: number;
    breakdown: MassBreakdownEntry[];
    /** True if any contributing element uses a most-stable-isotope mass. */
    usedMassNumber: boolean;
}
/** Sum a parsed formula into molar mass + per-element breakdown. */
export declare function formulaMass(parse: Extract<FormulaParse, {
    ok: true;
}>): MassResult;
/** Resolve a molar mass from either a parsed formula or an explicit number. */
export declare function resolveMolarMass(formula: string | undefined, molarMass: number | undefined): {
    mass: number;
    formula: string | null;
    error: string | null;
};
export interface SolutionPercentInput {
    percentType: 'wv' | 'vv' | 'mass';
    soluteG?: number;
    solutionMl?: number;
    soluteMl?: number;
    totalMl?: number;
    soluteMassG?: number;
    totalMassG?: number;
}
export type SolutionPercentResult = {
    ok: true;
    percent: number;
} | {
    ok: false;
    error: string;
};
/** Percent-solution math: w/v (g per 100 mL), v/v (mL per 100 mL) and mass %. */
export declare function solutionPercent(input: SolutionPercentInput): SolutionPercentResult;
export interface PhResult {
    pH: number;
    pOH: number;
    hConc: number;
    ohConc: number;
}
/** pH <-> concentration conversions at 25 °C (Kw = 1e-14). */
export declare function phFromH(hConc: number): PhResult;
export declare function phFromOH(ohConc: number): PhResult;
export declare function phFromPH(pH: number): PhResult;
export declare function phFromPOH(pOH: number): PhResult;
/** Henderson-Hasselbalch: pH = pKa + log10([A-]/[HA]). */
export declare function bufferPh(pKa: number, baseConc: number, acidConc: number): {
    pH: number;
    ratio: number;
};
//# sourceMappingURL=core.d.ts.map