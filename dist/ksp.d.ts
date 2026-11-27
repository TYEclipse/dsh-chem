/**
 * Solubility product (Ksp) → molar solubility, with the ion grouping made
 * explicit and the common-ion effect solved rather than approximated.
 *
 * The arithmetic of a 1:1 salt is trivial (`s = √Ksp`), and that is exactly why
 * the tool exists: the hard part is *which ions the formula actually breaks
 * into*. `Ca(OH)2` is one Ca²⁺ and two OH⁻ — not three "element counts" — and
 * `Ca3(PO4)2` gives `Ksp = (3s)³(2s)² = 108·s⁵`, a different exponent than the
 * 1:1 case. Getting the grouping wrong silently changes both the exponent and
 * the constant in front of it.
 *
 * Grouping rules (documented in the README):
 *   - the formula is split at one point into a cation half and an anion half,
 *     and each half must match a known ion exactly (element composition and
 *     charge) — a split that does not match is simply not a candidate;
 *   - a candidate must be electroneutral (`a·q_cation + b·q_anion = 0`) and its
 *     coefficients must be whole numbers;
 *   - **more than one candidate is reported as ambiguous, not resolved by
 *     preference** (`Hg2Cl2` is both `Hg2²⁺ + 2 Cl⁻` and `2 Hg⁺ + 2 Cl⁻`, which
 *     give different answers). Pass `cation` / `anion` to choose.
 *
 * Honest contracts: hydrates are rejected (Ksp describes the anhydrous solid),
 * a charged or coefficient-prefixed formula is rejected, an ion that is not in
 * the table is named in the error, and Ksp is used exactly as supplied — this
 * module does not interpolate temperature, activities or ionic strength.
 *
 * @module dsh-chem/ksp
 */
export interface KspIonOut {
    /** Ion symbol as printed, e.g. "Ca^2+". */
    symbol: string;
    /** Net charge in elementary units. */
    charge: number;
    /** How many of this ion one formula unit releases. */
    count: number;
    /** Concentration at saturation, mol/L (adds any common ion already present). */
    concentration: number;
}
export interface CommonIonOut {
    /** Ion symbol supplied as the common ion. */
    symbol: string;
    /** Concentration of that ion already present before dissolution, mol/L. */
    concentration: number;
    /** Molar solubility in the presence of that common ion. */
    molarSolubility: number;
    /** Mass solubility in the presence of that common ion, g/L. */
    massSolubility: number;
    /** Plain-water solubility ÷ common-ion solubility (≥ 1). */
    reductionFactor: number;
    /** Ksp recomputed from the solved concentrations (relative residual). */
    residual: number;
}
export interface KspSuccess {
    ok: true;
    formula: string;
    ksp: number;
    pKsp: number;
    cation: KspIonOut;
    anion: KspIonOut;
    /** a^a · b^b — the constant in front of s^(a+b). */
    stoichiometricFactor: number;
    /** a + b — the exponent on s. */
    ionCount: number;
    /** Printed form, e.g. "Ksp = [Ca^2+]^3·[PO4^3-]^2 = (3s)^3(2s)^2 = 108·s^5". */
    expression: string;
    molarMass: number;
    /** Plain-water molar solubility, mol/L. */
    molarSolubility: number;
    /** Plain-water mass solubility, g/L. */
    massSolubility: number;
    /** Charge balance check: a·q_cation + b·q_anion (must be 0). */
    chargeBalance: number;
    commonIon?: CommonIonOut;
    notes: string[];
}
export interface KspFailure {
    ok: false;
    error: string;
    /** True when the formula matches more than one ion grouping. */
    ambiguous?: boolean;
    candidates?: string[];
}
export type KspResult = KspSuccess | KspFailure;
export interface KspInput {
    formula: string;
    ksp: number;
    /** Force the cation, e.g. "Ca2+", "NH4+", "Hg2^2+". */
    cation?: string;
    /** Force the anion, e.g. "OH-", "PO4^3-". */
    anion?: string;
    /** Ion symbol (either of the salt's own ions) already present at this concentration. */
    commonIon?: string;
    /** Concentration of the common ion in mol/L (required with `commonIon`). */
    commonIonConcentration?: number;
}
/**
 * Ksp → molar solubility. See the module header for the grouping rules.
 */
export declare function kspSolubility(input: KspInput): KspResult;
//# sourceMappingURL=ksp.d.ts.map