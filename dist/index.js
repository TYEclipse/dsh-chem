/**
 * dsh-chem — deterministic chemistry math for DeepSeek Harness.
 *
 * Seven zero-dependency tools (pure arithmetic over IUPAC data):
 *   formula_mass      — parse a chemical formula and compute its molar mass + breakdown
 *   mass_moles        — convert mass (g) <-> amount (mol) <-> particle count
 *   solution          — molarity, dilution (C1V1 = C2V2) and percent solutions
 *   ph                — pH / pOH / [H+] / [OH-] conversions + Henderson-Hasselbalch buffers
 *   balance_equation  — balance an equation with exact rational arithmetic (+ charge balance)
 *   gas_law           — PV = nRT solves, molar volume, gas density and molar mass from density
 *   stoichiometry     — limiting reagent, extent of reaction, product mass and leftovers
 *
 * Agents get chemistry arithmetic wrong in predictable ways: molar-mass sums,
 * case-sensitivity of element symbols (Co vs CO), dilution factors, log/antilog
 * pH math and equation coefficients. These tools replace all of it with
 * deterministic tables and formulas.
 *
 * @module dsh-chem
 */
import z from '@deepseek-ai/schemastery';
import { buildChemTools } from "./tools.js";
/** Stable Cordis plugin name (also the config key under `plugins:`). */
export const name = 'dsh-chem';
/** Services required before tool registration can start. */
export const inject = ['agents', 'tools'];
export const Config = z.object({});
/** Register every chemistry tool on one agent; returns the disposer. */
function decorate(agent, tools) {
    const disposers = Object.values(tools).map((definition) => agent.ctx.tools.register(definition));
    return () => {
        for (const dispose of disposers) {
            try {
                dispose();
            }
            catch {
                // already disposed
            }
        }
    };
}
/** Mount the chemistry tools on every live agent and every future one. */
export function apply(ctx, _config) {
    const tools = buildChemTools();
    const disposers = new Set();
    const decorateAgent = (agent) => {
        try {
            disposers.add(decorate(agent, tools));
        }
        catch (error) {
            ctx.logger('chem').warn(`tool registration for agent ${agent.id} failed: ${error instanceof Error ? error.message : String(error)}`);
        }
    };
    for (const agent of ctx.agents.list())
        decorateAgent(agent);
    const off = ctx.on('agent/created', ({ agent }) => decorateAgent(agent));
    ctx.effect(() => () => {
        off();
        for (const dispose of disposers) {
            try {
                dispose();
            }
            catch {
                // already disposed
            }
        }
        disposers.clear();
    });
}
//# sourceMappingURL=index.js.map