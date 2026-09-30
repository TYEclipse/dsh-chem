/**
 * Stoichiometry core for dsh-chem: balance a chemical equation with exact
 * rational arithmetic (no floating point anywhere in the solver).
 *
 * Why exact: balancing is an integer nullspace problem. Doing it in doubles
 * produces coefficients like 2.9999999997 and "smallest integer solution"
 * claims that are simply wrong. Here every operation runs over BigInt-backed
 * fractions, then the solution is cleared of denominators, divided by the GCD
 * and sign-normalized.
 *
 * Deterministic rules (documented in the README):
 *   - arrow: exactly one of `->`, `→`, `=>`, `==>`, `=`; species on each side
 *   - species separator: a `+` that is surrounded by spaces, or a compact `+`
 *     that follows a letter/`)`/digit and is followed by an uppercase letter or
 *     `(` (so `Fe+O2` splits while `Na+ + Cl-` keeps its charges — write charged
 *     species in the spaced form, the compact form is ambiguous by construction)
 *   - a leading coefficient in a species (`2H2O`) is parsed and ignored: scaling
 *     one column cannot change the solution set, and the output always prints
 *     the smallest integer coefficients
 *   - at least two species are required; every species must end up with a
 *     positive integer coefficient, otherwise the equation does not balance
 *
 * @module dsh-chem/balance
 */
/** One formula unit on either side of the arrow. */
export interface BalanceSpecies {
    /** Normalized formula without a leading coefficient, e.g. `Fe2O3`. */
    formula: string;
    /** Element counts per formula unit, hydrate water included. */
    counts: Map<string, number>;
    /** Signed ionic charge parsed from the species (0 when neutral). */
    charge: number;
}
/** A balanced term: formula + smallest integer coefficient. */
export interface BalanceTerm {
    formula: string;
    coefficient: number;
}
/** Per-element atom conservation report. */
export interface ElementBalance {
    symbol: string;
    left: number;
    right: number;
}
/** Machine-readable failure kind, so callers can explain *why* nothing balances. */
export type BalanceFailureKind = 'parse' | 'inconsistent' | 'degenerate' | 'no-positive' | 'too-ambitious';
export type BalanceResult = {
    ok: true;
    reactants: BalanceTerm[];
    products: BalanceTerm[];
    elements: ElementBalance[];
    /** False when several independent balances exist (e.g. two reactions at once). */
    unique: boolean;
    /** True when at least one species carries a charge (charge row included). */
    chargeBalanced: boolean;
} | {
    ok: false;
    kind: BalanceFailureKind;
    error: string;
};
/** Balance `equation`, e.g. `"C3H8 + O2 -> CO2 + H2O"`. */
export declare function balanceEquation(equation: string): BalanceResult;
/** Render a balanced result as `2 Fe + 3 O2 -> Fe2O3`-style text. */
export declare function renderBalancedEquation(reactants: BalanceTerm[], products: BalanceTerm[]): string;
//# sourceMappingURL=balance.d.ts.map