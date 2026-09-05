/**
 * Pure chemistry math for dsh-chem: formula parsing, molar mass, mole/particle
 * conversions, solution concentration and pH — zero dependencies, no I/O.
 *
 * All rules are deterministic and documented; agents should call these tools
 * instead of doing molar-mass / dilution / pH arithmetic in their head.
 *
 * @module dsh-chem/core
 */
import { elementBySymbolLenient } from "./periodic.js";
/** SI Avogadro constant (exact, mol^-1, CODATA 2018). */
export const AVOGADRO = 6.02214076e23;
/** Water autoionization constant at 25 °C (mol²/L²). */
export const KW_25C = 1.0e-14;
/** Hydrate separators accepted by the parser (middle dot and asterisk). */
const HYDRATE_SEPARATORS = new Set(['·', '*', '•']);
function isDigit(ch) {
    return ch >= '0' && ch <= '9';
}
/** Read a run of digits at the current position; returns the integer (0 if none). */
function readInteger(state) {
    let value = 0;
    while (state.pos < state.input.length && isDigit(state.input[state.pos])) {
        value = value * 10 + (state.input.charCodeAt(state.pos) - 48);
        state.pos += 1;
    }
    return value;
}
function skipWhitespace(state) {
    while (state.pos < state.input.length && /\s/.test(state.input[state.pos])) {
        state.pos += 1;
    }
}
/** Parse one atom or parenthesized group (without multiplicity); returns null on error. */
function parseUnit(state) {
    skipWhitespace(state);
    const ch = state.input[state.pos];
    if (ch === undefined)
        return null;
    if (ch === '(') {
        state.pos += 1;
        const inner = [];
        let innerSource = '(';
        for (;;) {
            const unit = parseUnit(state);
            if (unit === null)
                return null;
            // Inside-paren multiplicity, e.g. (CH3)2CO -> the 3 is H's subscript.
            if (isDigit(state.input[state.pos])) {
                const innerMult = readInteger(state);
                if (innerMult === 0)
                    return null;
                inner.push(...unit.atoms.map((a) => ({ symbol: a.symbol, count: a.count * innerMult })));
                innerSource += unit.source + String(innerMult);
            }
            else {
                inner.push(...unit.atoms);
                innerSource += unit.source;
            }
            skipWhitespace(state);
            const next = state.input[state.pos];
            if (next === ')') {
                state.pos += 1;
                innerSource += ')';
                return { atoms: inner, source: innerSource };
            }
            if (next === undefined)
                return null;
            // Not a closing paren: must be the start of the next unit (loop again).
        }
    }
    if (ch >= 'A' && ch <= 'Z') {
        state.pos += 1;
        let symbol = ch;
        const lower = state.input[state.pos];
        if (lower !== undefined && lower >= 'a' && lower <= 'z') {
            symbol += lower;
            state.pos += 1;
        }
        return { atoms: [{ symbol, count: 1 }], source: symbol };
    }
    if (ch >= 'a' && ch <= 'z') {
        // Lenient: a single lowercase letter that is a real element symbol is
        // accepted and normalized ("h2o" -> H2O, "nacl" fails on 'a' with guidance).
        // Multi-letter lowercase runs are NOT guessed ("co" is ambiguous: cobalt
        // Co vs carbon+oxygen C O) — those fail with the case-sensitivity hint.
        const candidate = elementBySymbolLenient(ch);
        if (candidate && candidate.symbol.toLowerCase() === ch) {
            state.pos += 1;
            return { atoms: [{ symbol: candidate.symbol, count: 1 }], source: candidate.symbol };
        }
        return null;
    }
    return null;
}
/**
 * Parse a chemical formula: nested parentheses, hydrate dot (CuSO4·5H2O),
 * leading coefficient (2H2O) and trailing ionic charge (Fe3+, SO4^2-, OH-).
 */
export function parseFormula(input) {
    const raw = input.trim();
    if (raw.length === 0)
        return { ok: false, error: 'formula is empty' };
    if (raw.length > 400)
        return { ok: false, error: `formula is too long (${raw.length} chars; max 400)` };
    const state = { input: raw, pos: 0 };
    // Leading coefficient.
    const coefficient = isDigit(state.input[0]) ? readInteger(state) : null;
    if (coefficient !== null && coefficient === 0) {
        return { ok: false, error: 'leading coefficient must be positive (got 0)' };
    }
    if (coefficient !== null && state.pos === raw.length) {
        return { ok: false, error: `formula "${raw}" has a coefficient but no formula` };
    }
    // Main formula groups.
    const groups = [];
    let hydrateGroups = null;
    let charge = null;
    const parseSequence = () => {
        for (;;) {
            skipWhitespace(state);
            const ch = state.input[state.pos];
            if (ch === undefined)
                return true;
            if (HYDRATE_SEPARATORS.has(ch) || ch === '+' || ch === '-' || ch === '^')
                return true;
            const unit = parseUnit(state);
            if (unit === null) {
                return false;
            }
            const digitStart = state.pos;
            const multiplicity = isDigit(state.input[state.pos]) ? readInteger(state) : 1;
            if (multiplicity === 0) {
                return false;
            }
            // "Fe3+" vs "NH4+" ambiguity: digits followed by a +/- sign belong to the
            // IONIC CHARGE only when every element parsed so far is the same single
            // element ("Fe3+" -> Fe with 3+, "S2-" -> S with 2-, "Ca2+" -> Ca with
            // 2+). Otherwise they are a subscript and the bare sign is a 1-unit
            // charge ("NH4+" -> N + H4 with +, "NO3-" -> N + O3 with -). The
            // explicit "^" form is always the charge ("SO4^2-", "O2^+").
            if (state.pos > digitStart) {
                const after = state.input[state.pos];
                if (after === '+' || after === '-') {
                    const distinct = new Set();
                    for (const g of groups)
                        for (const a of g.atoms)
                            distinct.add(a.symbol);
                    for (const a of unit.atoms)
                        distinct.add(a.symbol);
                    if (distinct.size === 1) {
                        // Digits belong to the charge: keep the unit at count 1, rewind so
                        // the charge parser consumes digits+sign ("Fe3+" -> Fe with 3+).
                        groups.push({ atoms: unit.atoms, source: unit.source });
                        state.pos = digitStart;
                        return true;
                    }
                }
            }
            const atoms = unit.atoms.map((a) => ({ symbol: a.symbol, count: a.count * multiplicity }));
            groups.push({ atoms, source: unit.source + (multiplicity > 1 ? String(multiplicity) : '') });
            if (groups.length > 400)
                return false;
        }
    };
    if (!parseSequence()) {
        return { ok: false, error: `could not parse formula "${raw}" near position ${state.pos}` };
    }
    // Hydrate part after a separator.
    skipWhitespace(state);
    const sep = state.input[state.pos];
    let hydrateCoefficient = 1;
    if (sep !== undefined && HYDRATE_SEPARATORS.has(sep)) {
        state.pos += 1;
        skipWhitespace(state);
        // Hydrate coefficient (water molecules per formula unit), e.g. 5 in ·5H2O.
        hydrateCoefficient = 1;
        if (isDigit(state.input[state.pos])) {
            hydrateCoefficient = readInteger(state);
            if (hydrateCoefficient === 0) {
                return { ok: false, error: `hydrate coefficient must be positive in "${raw}"` };
            }
        }
        hydrateGroups = [];
        for (;;) {
            skipWhitespace(state);
            const ch = state.input[state.pos];
            if (ch === undefined)
                break;
            if (ch === '+' || ch === '-' || ch === '^' || HYDRATE_SEPARATORS.has(ch))
                break;
            const unit = parseUnit(state);
            if (unit === null) {
                return { ok: false, error: `could not parse hydrate part near position ${state.pos}` };
            }
            const multiplicity = isDigit(state.input[state.pos]) ? readInteger(state) : 1;
            if (multiplicity === 0) {
                return { ok: false, error: `hydrate multiplicity must be positive in "${raw}"` };
            }
            hydrateGroups.push({
                atoms: unit.atoms.map((a) => ({ symbol: a.symbol, count: a.count * multiplicity * hydrateCoefficient })),
                source: unit.source + (multiplicity > 1 ? String(multiplicity) : ''),
            });
        }
    }
    // Trailing ionic charge: optional '^', digits, then sign.
    skipWhitespace(state);
    let chargeStr = null;
    {
        let pos = state.pos;
        let magnitude = '';
        if (state.input[pos] === '^')
            pos += 1;
        while (pos < raw.length && isDigit(raw[pos])) {
            magnitude += raw[pos];
            pos += 1;
        }
        if (pos < raw.length && (raw[pos] === '+' || raw[pos] === '-')) {
            chargeStr = magnitude + raw[pos];
            state.pos = pos + 1;
        }
    }
    if (chargeStr !== null) {
        charge = chargeStr;
        skipWhitespace(state);
        if (state.pos !== raw.length) {
            return { ok: false, error: `unexpected trailing content in "${raw}" after charge "${chargeStr}"` };
        }
    }
    else if (state.pos !== raw.length) {
        return { ok: false, error: `unexpected trailing content in "${raw}" near position ${state.pos}` };
    }
    // Validate every symbol against the table (case-sensitive first, lenient fallback).
    const allGroups = hydrateGroups !== null ? [...groups, ...hydrateGroups] : groups;
    for (const group of allGroups) {
        for (const atom of group.atoms) {
            if (!elementBySymbolLenient(atom.symbol)) {
                return { ok: false, error: `unknown element symbol "${atom.symbol}" in "${raw}" (symbols are case-sensitive: "Co" is cobalt, "CO" is carbon+oxygen)` };
            }
        }
    }
    // totals covers the MAIN groups only; hydrate is reported separately and
    // added by formulaMass (water molecules per formula unit).
    const totals = mergeGroups(groups);
    const hydrate = hydrateGroups === null ? null : mergeGroups(hydrateGroups);
    const normalized = (coefficient !== null ? String(coefficient) : '')
        + groups.map((g) => g.source).join('')
        + (hydrateGroups !== null
            ? '·' + (hydrateCoefficient > 1 ? String(hydrateCoefficient) : '')
                + (hydrateGroups.length > 0 ? hydrateGroups.map((g) => g.source).join('') : '')
            : '')
        + (charge ?? '');
    return {
        ok: true,
        groups: groups.map((g) => ({ elements: mergeGroups([g]), source: g.source })),
        totals,
        hydrate,
        coefficient,
        charge,
        normalized,
    };
}
function mergeGroups(groups) {
    const seen = new Map();
    for (const group of groups) {
        for (const atom of group.atoms) {
            seen.set(atom.symbol, (seen.get(atom.symbol) ?? 0) + atom.count);
        }
    }
    return [...seen.entries()].map(([symbol, count]) => ({ symbol, count }));
}
/** Sum a parsed formula into molar mass + per-element breakdown. */
export function formulaMass(parse) {
    const contrib = new Map();
    const add = (counts) => {
        for (const c of counts) {
            const elem = elementBySymbolLenient(c.symbol);
            const prev = contrib.get(c.symbol);
            if (prev)
                prev.count += c.count;
            else
                contrib.set(c.symbol, { elem, count: c.count });
        }
    };
    add(parse.totals);
    if (parse.hydrate)
        add(parse.hydrate);
    let total = 0;
    const rows = [];
    for (const { elem, count } of contrib.values()) {
        const mass = elem.atomicMass * count;
        total += mass;
        rows.push({ symbol: elem.symbol, count, elem, mass });
    }
    rows.sort((a, b) => a.elem.atomicNumber - b.elem.atomicNumber);
    return {
        molarMass: total,
        usedMassNumber: rows.some((r) => r.elem.massNumber),
        breakdown: rows.map((r) => ({
            symbol: r.symbol,
            name: r.elem.name,
            count: r.count,
            atomicMass: r.elem.atomicMass,
            massContribution: r.mass,
            massFraction: total > 0 ? r.mass / total : 0,
        })),
    };
}
/** Resolve a molar mass from either a parsed formula or an explicit number. */
export function resolveMolarMass(formula, molarMass) {
    if (formula !== undefined && formula.trim() !== '') {
        const parsed = parseFormula(formula);
        if (!parsed.ok)
            return { mass: NaN, formula: null, error: parsed.error };
        return { mass: formulaMass(parsed).molarMass, formula: parsed.normalized, error: null };
    }
    if (molarMass === undefined) {
        return { mass: NaN, formula: null, error: 'provide either formula or molarMass' };
    }
    if (!Number.isFinite(molarMass) || molarMass <= 0) {
        return { mass: NaN, formula: null, error: `molarMass must be a positive number (got ${molarMass})` };
    }
    return { mass: molarMass, formula: null, error: null };
}
/** Percent-solution math: w/v (g per 100 mL), v/v (mL per 100 mL) and mass %. */
export function solutionPercent(input) {
    if (input.percentType === 'wv') {
        if (input.soluteG === undefined || input.solutionMl === undefined) {
            return { ok: false, error: 'percentType "wv" requires soluteG (grams) and solutionMl (milliliters)' };
        }
        if (!Number.isFinite(input.soluteG) || input.soluteG < 0)
            return { ok: false, error: `soluteG must be >= 0 (got ${input.soluteG})` };
        if (!Number.isFinite(input.solutionMl) || input.solutionMl <= 0)
            return { ok: false, error: `solutionMl must be positive (got ${input.solutionMl})` };
        return { ok: true, percent: (input.soluteG / input.solutionMl) * 100 };
    }
    if (input.percentType === 'vv') {
        if (input.soluteMl === undefined || input.totalMl === undefined) {
            return { ok: false, error: 'percentType "vv" requires soluteMl and totalMl (milliliters)' };
        }
        if (!Number.isFinite(input.soluteMl) || input.soluteMl < 0)
            return { ok: false, error: `soluteMl must be >= 0 (got ${input.soluteMl})` };
        if (!Number.isFinite(input.totalMl) || input.totalMl <= 0)
            return { ok: false, error: `totalMl must be positive (got ${input.totalMl})` };
        return { ok: true, percent: (input.soluteMl / input.totalMl) * 100 };
    }
    if (input.soluteMassG === undefined || input.totalMassG === undefined) {
        return { ok: false, error: 'percentType "mass" requires soluteMassG and totalMassG (grams)' };
    }
    if (!Number.isFinite(input.soluteMassG) || input.soluteMassG < 0)
        return { ok: false, error: `soluteMassG must be >= 0 (got ${input.soluteMassG})` };
    if (!Number.isFinite(input.totalMassG) || input.totalMassG <= 0)
        return { ok: false, error: `totalMassG must be positive (got ${input.totalMassG})` };
    return { ok: true, percent: (input.soluteMassG / input.totalMassG) * 100 };
}
/** pH <-> concentration conversions at 25 °C (Kw = 1e-14). */
export function phFromH(hConc) {
    const pH = -Math.log10(hConc);
    const pOH = 14 - pH;
    return { pH, pOH, hConc, ohConc: Math.pow(10, -pOH) };
}
export function phFromOH(ohConc) {
    const pOH = -Math.log10(ohConc);
    const pH = 14 - pOH;
    return { pH, pOH, hConc: Math.pow(10, -pH), ohConc };
}
export function phFromPH(pH) {
    const hConc = Math.pow(10, -pH);
    const pOH = 14 - pH;
    return { pH, pOH, hConc, ohConc: Math.pow(10, -pOH) };
}
export function phFromPOH(pOH) {
    const ohConc = Math.pow(10, -pOH);
    const pH = 14 - pOH;
    return { pH, pOH, hConc: Math.pow(10, -pH), ohConc };
}
/** Henderson-Hasselbalch: pH = pKa + log10([A-]/[HA]). */
export function bufferPh(pKa, baseConc, acidConc) {
    const ratio = baseConc / acidConc;
    return { pH: pKa + Math.log10(ratio), ratio };
}
//# sourceMappingURL=core.js.map