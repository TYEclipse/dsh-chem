/**
 * Stoichiometry on top of the exact-coefficient balancer: given how much of
 * each reactant you actually have, which one runs out first, how much product
 * that yields, and what is left over.
 *
 * Every number here is either a user input or derived from the balanced
 * coefficients and the IUPAC molar masses — nothing is estimated.
 *
 * @module dsh-chem/stoich
 */
export interface StoichInput {
    equation: string;
    /** One entry per reactant in the equation; give exactly one of moles / grams. */
    reactants: {
        formula: string;
        moles?: number;
        grams?: number;
    }[];
    /** Product whose theoretical yield should be reported (required for percent yield). */
    product?: string;
    /** Measured mass of `product`, in grams — turns the answer into a percent yield. */
    actualYieldG?: number;
}
export interface StoichReactantOut {
    formula: string;
    coefficient: number;
    molarMass: number;
    availableMoles: number;
    availableGrams: number;
    consumedMoles: number;
    leftoverMoles: number;
    leftoverGrams: number;
    limiting: boolean;
}
export interface StoichProductOut {
    formula: string;
    coefficient: number;
    molarMass: number;
    moles: number;
    grams: number;
}
export type StoichResult = {
    ok: true;
    equation: string;
    extent: number;
    limiting: string[];
    reactants: StoichReactantOut[];
    products: StoichProductOut[];
    yieldProduct?: string;
    theoreticalYieldG?: number;
    theoreticalYieldMol?: number;
    actualYieldG?: number;
    percentYield?: number;
} | {
    ok: false;
    error: string;
};
export declare function stoichiometry(input: StoichInput): StoichResult;
//# sourceMappingURL=stoich.d.ts.map