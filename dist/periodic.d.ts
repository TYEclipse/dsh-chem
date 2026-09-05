/**
 * Periodic table data for dsh-chem.
 *
 * Atomic masses are the IUPAC standard atomic weights (abridged, conventional
 * values, widely used in textbooks). Elements with no stable isotopes use the
 * mass number of the most stable known isotope (marked with `massNumber`),
 * matching the bracketed-value convention of IUPAC tables.
 *
 * @module dsh-chem/periodic
 */
export interface ElementData {
    /** Element symbol, case-sensitive ("Co" is cobalt, "CO" is carbon monoxide). */
    symbol: string;
    name: string;
    /** Atomic number Z (1–118). */
    atomicNumber: number;
    /** Standard atomic weight in g/mol. */
    atomicMass: number;
    /** True when the mass is the most-stable-isotope mass number (radioelement). */
    massNumber: boolean;
}
/** All 118 elements, ordered by atomic number. */
export declare const ELEMENTS: readonly ElementData[];
export declare function elementBySymbol(symbol: string): ElementData | undefined;
/**
 * Lenient lookup: exact symbol first, then lowercase-insensitive match
 * ("co" -> Co cobalt). This catches the most common LLM mistake (lowercasing
 * symbols) while still treating "CO" as carbon + oxygen.
 */
export declare function elementBySymbolLenient(symbol: string): ElementData | undefined;
//# sourceMappingURL=periodic.d.ts.map