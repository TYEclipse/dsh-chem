/**
 * dsh-chem — deterministic chemistry math for DeepSeek Harness.
 *
 * Four zero-dependency tools (pure arithmetic over IUPAC data):
 *   formula_mass — parse a chemical formula and compute its molar mass + breakdown
 *   mass_moles   — convert mass (g) <-> amount (mol) <-> particle count
 *   solution     — molarity, dilution (C1V1 = C2V2) and percent solutions
 *   ph           — pH / pOH / [H+] / [OH-] conversions + Henderson-Hasselbalch buffers
 *
 * Agents get chemistry arithmetic wrong in predictable ways: molar-mass sums,
 * case-sensitivity of element symbols (Co vs CO), dilution factors and log/antilog
 * pH math. These tools replace all of it with deterministic tables and formulas.
 *
 * @module dsh-chem
 */
import type { Context } from '@deepseek-ai/cordis';
import z from '@deepseek-ai/schemastery';
/** Stable Cordis plugin name (also the config key under `plugins:`). */
export declare const name = "dsh-chem";
/** Services required before tool registration can start. */
export declare const inject: string[];
/** Plugin configuration (no tunables today — kept for schema completeness). */
export interface Config {
    /** Reserved for future options. */
    _?: never;
}
export declare const Config: z<Config>;
/** Mount the chemistry tools on every live agent and every future one. */
export declare function apply(ctx: Context, _config: Config): void;
//# sourceMappingURL=index.d.ts.map