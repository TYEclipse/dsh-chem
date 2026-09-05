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
export declare function buildChemTools(): ToolSet;
//# sourceMappingURL=tools.d.ts.map