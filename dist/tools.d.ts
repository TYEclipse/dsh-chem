/**
 * Tool definitions for dsh-chem: four deterministic chemistry tools exposed to
 * every agent via defineTool. Strict JSON-schema parameter surfaces, explicit
 * result interfaces matching the inferred output types.
 *
 * @module dsh-chem/tools
 */
import { type ToolDefinition } from '@deepseek-ai/dsh-tools';
export interface ToolSet {
    formula_mass: ToolDefinition;
    mass_moles: ToolDefinition;
    solution: ToolDefinition;
    ph: ToolDefinition;
    balance_equation: ToolDefinition;
    gas_law: ToolDefinition;
    stoichiometry: ToolDefinition;
}
/** Output-facing shape of `stoichiometry` (undefined keys are never emitted). */
export interface StoichiometryOut {
    valid: boolean;
    equation?: string;
    extent?: number;
    limiting?: string[];
    reactants?: {
        formula: string;
        coefficient: number;
        molarMass: number;
        availableMoles: number;
        availableGrams: number;
        consumedMoles: number;
        leftoverMoles: number;
        leftoverGrams: number;
        limiting: boolean;
    }[];
    products?: {
        formula: string;
        coefficient: number;
        molarMass: number;
        moles: number;
        grams: number;
    }[];
    yieldProduct?: string;
    theoreticalYieldMol?: number;
    theoreticalYieldG?: number;
    actualYieldG?: number;
    percentYield?: number;
    error?: string;
}
export interface FormulaMassResult {
    valid: boolean;
    formula?: string;
    molarMass?: number;
    elements?: {
        symbol: string;
        name: string;
        count: number;
        atomicMass: number;
        massContribution: number;
        massFraction: number;
    }[];
    hydrate?: {
        symbol: string;
        count: number;
    }[];
    charge?: string;
    coefficient?: number;
    usedMassNumber?: boolean;
    error?: string;
}
export interface MassMolesResult {
    valid: boolean;
    mode?: string;
    moles?: number;
    massG?: number;
    particles?: number;
    molarMass?: number;
    formula?: string;
    avogadro?: number;
    error?: string;
}
export interface SolutionResult {
    valid: boolean;
    operation?: string;
    molarity?: number;
    massG?: number;
    volumeMl?: number;
    molarMass?: number;
    formula?: string;
    c1?: number;
    v1?: number;
    c2?: number;
    v2?: number;
    solved?: string;
    percent?: number;
    percentType?: string;
    error?: string;
}
export interface PhResultOut {
    valid: boolean;
    operation?: string;
    pH?: number;
    pOH?: number;
    hConc?: number;
    ohConc?: number;
    bufferRatio?: number;
    pKa?: number;
    error?: string;
}
export interface BalanceTermOut {
    formula: string;
    coefficient: number;
}
export interface BalanceEquationResult {
    valid: boolean;
    /** false when the species cannot balance at all (valid stays true: the input was understood). */
    balanced?: boolean;
    /** false when several independent balances exist; the smallest one is returned. */
    unique?: boolean;
    equation?: string;
    reactants?: BalanceTermOut[];
    products?: BalanceTermOut[];
    elements?: {
        symbol: string;
        left: number;
        right: number;
    }[];
    /** true when at least one species carried a charge (charge conservation checked too). */
    chargeBalanced?: boolean;
    error?: string;
}
export interface GasLawResult {
    valid: boolean;
    operation?: string;
    solved?: string;
    pressureAtm?: number;
    pressureKpa?: number;
    pressureMmhg?: number;
    volumeL?: number;
    volumeMl?: number;
    moles?: number;
    temperatureK?: number;
    temperatureC?: number;
    molarVolumeLPerMol?: number;
    molarVolumeMlPerMol?: number;
    densityGPerL?: number;
    densityGPerMl?: number;
    molarMass?: number;
    formula?: string;
    r?: number;
    error?: string;
}
export declare function buildChemTools(): ToolSet;
//# sourceMappingURL=tools.d.ts.map