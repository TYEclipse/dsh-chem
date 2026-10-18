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
import { balanceEquation } from "./balance.js";
import { formulaMass, parseFormula } from "./core.js";
/** Fold `-0` and binary-noise tails at the output boundary (never mid-computation). */
function roundOut(value, digits = 12) {
    if (!Number.isFinite(value))
        return value;
    const r = Number(value.toPrecision(digits));
    return r === 0 ? 0 : r;
}
function molarMassOf(formula) {
    const parsed = parseFormula(formula);
    if (!parsed.ok)
        return { error: parsed.error };
    return formulaMass(parsed).molarMass;
}
export function stoichiometry(input) {
    const balanced = balanceEquation(input.equation);
    if (!balanced.ok)
        return { ok: false, error: balanced.error };
    const equationText = [
        balanced.reactants.map((t) => `${t.coefficient === 1 ? '' : t.coefficient}${t.formula}`).join(' + '),
        balanced.products.map((t) => `${t.coefficient === 1 ? '' : t.coefficient}${t.formula}`).join(' + '),
    ].join(' -> ');
    // ---- normalize + validate the supplied amounts -------------------------
    const supplied = new Map();
    const given = input.reactants ?? [];
    if (given.length === 0)
        return { ok: false, error: 'reactants must list the amount of every reactant in the equation' };
    for (const entry of given) {
        const parsed = parseFormula(entry.formula);
        if (!parsed.ok)
            return { ok: false, error: `reactant "${entry.formula}": ${parsed.error}` };
        const formula = parsed.normalized;
        if (supplied.has(formula))
            return { ok: false, error: `reactant "${formula}" is listed more than once` };
        const hasMoles = entry.moles !== undefined;
        const hasGrams = entry.grams !== undefined;
        if (hasMoles === hasGrams) {
            return { ok: false, error: `reactant "${formula}" needs exactly one of moles or grams (got ${hasMoles ? 'both' : 'neither'})` };
        }
        let moles;
        if (hasMoles) {
            if (!Number.isFinite(entry.moles) || entry.moles <= 0) {
                return { ok: false, error: `reactant "${formula}": moles must be a positive number (got ${entry.moles})` };
            }
            moles = entry.moles;
        }
        else {
            if (!Number.isFinite(entry.grams) || entry.grams <= 0) {
                return { ok: false, error: `reactant "${formula}": grams must be a positive number (got ${entry.grams})` };
            }
            const mass = molarMassOf(formula);
            if (typeof mass !== 'number')
                return { ok: false, error: `reactant "${formula}": ${mass.error}` };
            moles = entry.grams / mass;
        }
        supplied.set(formula, { formula, moles });
    }
    const known = balanced.reactants.map((t) => t.formula);
    for (const formula of supplied.keys()) {
        if (!known.includes(formula)) {
            return { ok: false, error: `"${formula}" is not a reactant in this equation (reactants: ${known.join(', ')})` };
        }
    }
    const missing = known.filter((formula) => !supplied.has(formula));
    if (missing.length > 0) {
        return { ok: false, error: `no amount given for reactant(s): ${missing.join(', ')} — every reactant needs one, otherwise the limiting reagent is undefined` };
    }
    // ---- extent of reaction ------------------------------------------------
    const ratios = balanced.reactants.map((term) => supplied.get(term.formula).moles / term.coefficient);
    const extent = Math.min(...ratios);
    const tol = 1e-9 * Math.max(1, Math.abs(extent));
    const limiting = ratios
        .map((ratio, index) => (Math.abs(ratio - extent) <= tol ? balanced.reactants[index].formula : null))
        .filter((formula) => formula !== null);
    const reactants = balanced.reactants.map((term) => {
        const { moles } = supplied.get(term.formula);
        const mass = molarMassOf(term.formula);
        const molarMass = typeof mass === 'number' ? mass : Number.NaN;
        const consumed = extent * term.coefficient;
        const leftover = Math.max(0, moles - consumed);
        return {
            formula: term.formula,
            coefficient: term.coefficient,
            molarMass: roundOut(molarMass),
            availableMoles: roundOut(moles),
            availableGrams: roundOut(moles * molarMass),
            consumedMoles: roundOut(consumed),
            leftoverMoles: roundOut(leftover),
            leftoverGrams: roundOut(leftover * molarMass),
            limiting: limiting.includes(term.formula),
        };
    });
    const products = balanced.products.map((term) => {
        const mass = molarMassOf(term.formula);
        const molarMass = typeof mass === 'number' ? mass : Number.NaN;
        const moles = extent * term.coefficient;
        return {
            formula: term.formula,
            coefficient: term.coefficient,
            molarMass: roundOut(molarMass),
            moles: roundOut(moles),
            grams: roundOut(moles * molarMass),
        };
    });
    const result = {
        ok: true,
        equation: equationText,
        extent: roundOut(extent),
        limiting,
        reactants,
        products,
    };
    if (input.product !== undefined) {
        const parsed = parseFormula(input.product);
        if (!parsed.ok)
            return { ok: false, error: `product "${input.product}": ${parsed.error}` };
        const target = products.find((p) => p.formula === parsed.normalized);
        if (!target) {
            return { ok: false, error: `"${parsed.normalized}" is not a product in this equation (products: ${products.map((p) => p.formula).join(', ')})` };
        }
        result.yieldProduct = target.formula;
        result.theoreticalYieldMol = target.moles;
        result.theoreticalYieldG = target.grams;
        if (input.actualYieldG !== undefined) {
            if (!Number.isFinite(input.actualYieldG) || input.actualYieldG < 0) {
                return { ok: false, error: `actualYieldG must be a non-negative number (got ${input.actualYieldG})` };
            }
            if (target.grams <= 0) {
                return { ok: false, error: 'theoretical yield is zero — a percent yield has no meaning here' };
            }
            result.actualYieldG = roundOut(input.actualYieldG);
            result.percentYield = roundOut((input.actualYieldG / target.grams) * 100);
        }
    }
    else if (input.actualYieldG !== undefined) {
        return { ok: false, error: 'actualYieldG needs a product: name the product whose yield you measured' };
    }
    return result;
}
//# sourceMappingURL=stoich.js.map