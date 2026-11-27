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
import { parseFormula, roundOut } from "./core.js";
import { formulaMass } from "./core.js";
/** Cations recognised when splitting a salt formula. */
const CATIONS = [
    { label: 'Li^+', formula: 'Li', charge: 1 },
    { label: 'Na^+', formula: 'Na', charge: 1 },
    { label: 'K^+', formula: 'K', charge: 1 },
    { label: 'Rb^+', formula: 'Rb', charge: 1 },
    { label: 'Cs^+', formula: 'Cs', charge: 1 },
    { label: 'Ag^+', formula: 'Ag', charge: 1 },
    { label: 'Cu^+', formula: 'Cu', charge: 1 },
    { label: 'Hg^+', formula: 'Hg', charge: 1 },
    { label: 'NH4^+', formula: 'NH4', charge: 1 },
    { label: 'Be^2+', formula: 'Be', charge: 2 },
    { label: 'Mg^2+', formula: 'Mg', charge: 2 },
    { label: 'Ca^2+', formula: 'Ca', charge: 2 },
    { label: 'Sr^2+', formula: 'Sr', charge: 2 },
    { label: 'Ba^2+', formula: 'Ba', charge: 2 },
    { label: 'Zn^2+', formula: 'Zn', charge: 2 },
    { label: 'Cd^2+', formula: 'Cd', charge: 2 },
    { label: 'Hg^2+', formula: 'Hg', charge: 2 },
    { label: 'Hg2^2+', formula: 'Hg2', charge: 2 },
    { label: 'Cu^2+', formula: 'Cu', charge: 2 },
    { label: 'Fe^2+', formula: 'Fe', charge: 2 },
    { label: 'Co^2+', formula: 'Co', charge: 2 },
    { label: 'Ni^2+', formula: 'Ni', charge: 2 },
    { label: 'Mn^2+', formula: 'Mn', charge: 2 },
    { label: 'Pb^2+', formula: 'Pb', charge: 2 },
    { label: 'Sn^2+', formula: 'Sn', charge: 2 },
    { label: 'Al^3+', formula: 'Al', charge: 3 },
    { label: 'Ga^3+', formula: 'Ga', charge: 3 },
    { label: 'In^3+', formula: 'In', charge: 3 },
    { label: 'Cr^3+', formula: 'Cr', charge: 3 },
    { label: 'Fe^3+', formula: 'Fe', charge: 3 },
    { label: 'Bi^3+', formula: 'Bi', charge: 3 },
    { label: 'Sn^4+', formula: 'Sn', charge: 4 },
];
/** Anions recognised when splitting a salt formula. */
const ANIONS = [
    { label: 'F^-', formula: 'F', charge: -1 },
    { label: 'Cl^-', formula: 'Cl', charge: -1 },
    { label: 'Br^-', formula: 'Br', charge: -1 },
    { label: 'I^-', formula: 'I', charge: -1 },
    { label: 'OH^-', formula: 'OH', charge: -1 },
    { label: 'NO3^-', formula: 'NO3', charge: -1 },
    { label: 'NO2^-', formula: 'NO2', charge: -1 },
    { label: 'HCO3^-', formula: 'CHO3', charge: -1 },
    { label: 'HSO4^-', formula: 'HSO4', charge: -1 },
    { label: 'H2PO4^-', formula: 'H2PO4', charge: -1 },
    { label: 'ClO^-', formula: 'ClO', charge: -1 },
    { label: 'ClO3^-', formula: 'ClO3', charge: -1 },
    { label: 'ClO4^-', formula: 'ClO4', charge: -1 },
    { label: 'MnO4^-', formula: 'MnO4', charge: -1 },
    { label: 'CN^-', formula: 'CN', charge: -1 },
    { label: 'SCN^-', formula: 'CNS', charge: -1 },
    { label: 'CH3COO^-', formula: 'C2H3O2', charge: -1 },
    { label: 'IO3^-', formula: 'IO3', charge: -1 },
    { label: 'BrO3^-', formula: 'BrO3', charge: -1 },
    { label: 'O^2-', formula: 'O', charge: -2 },
    { label: 'S^2-', formula: 'S', charge: -2 },
    { label: 'Se^2-', formula: 'Se', charge: -2 },
    { label: 'O2^2-', formula: 'O2', charge: -2 },
    { label: 'SO4^2-', formula: 'SO4', charge: -2 },
    { label: 'SO3^2-', formula: 'SO3', charge: -2 },
    { label: 'S2O3^2-', formula: 'S2O3', charge: -2 },
    { label: 'CO3^2-', formula: 'CO3', charge: -2 },
    { label: 'HPO4^2-', formula: 'HPO4', charge: -2 },
    { label: 'C2O4^2-', formula: 'C2O4', charge: -2 },
    { label: 'CrO4^2-', formula: 'CrO4', charge: -2 },
    { label: 'Cr2O7^2-', formula: 'Cr2O7', charge: -2 },
    { label: 'SiO3^2-', formula: 'SiO3', charge: -2 },
    { label: 'PO4^3-', formula: 'PO4', charge: -3 },
    { label: 'PO3^3-', formula: 'PO3', charge: -3 },
    { label: 'AsO4^3-', formula: 'AsO4', charge: -3 },
];
/** Element composition of an ion's neutral formula, as a sorted symbol→count map. */
function compositionOf(formula) {
    const parse = parseFormula(formula);
    if (!parse.ok)
        return null;
    const map = new Map();
    for (const entry of parse.totals)
        map.set(entry.symbol, entry.count);
    return map;
}
const ION_TABLE = [...CATIONS, ...ANIONS].flatMap((ion) => {
    const composition = compositionOf(ion.formula);
    return composition === null ? [] : [{ ion, composition }];
});
function sameShape(a, b) {
    if (a.size !== b.size)
        return false;
    for (const key of a.keys())
        if (!b.has(key))
            return false;
    return true;
}
/**
 * If `half` is `k` copies of `ion`'s composition (k ≥ 1 integer), return k.
 * Otherwise null.
 */
function multipleOf(half, ion) {
    if (!sameShape(half, ion))
        return null;
    let k = null;
    for (const [symbol, count] of ion) {
        const have = half.get(symbol) ?? 0;
        if (count <= 0 || have % count !== 0)
            return null;
        const ratio = have / count;
        if (k === null)
            k = ratio;
        else if (k !== ratio)
            return null;
    }
    return k !== null && k >= 1 ? k : null;
}
/** Half-formula element totals (groups already carry their multiplicities). */
function halfTotals(parse, from, to) {
    const map = new Map();
    for (let i = from; i < to; i += 1) {
        for (const element of parse.groups[i].elements) {
            map.set(element.symbol, (map.get(element.symbol) ?? 0) + element.count);
        }
    }
    return map;
}
/** Format a candidate grouping for the ambiguity message. */
function describe(cation, anion) {
    return `${cation.count} ${cation.ion.label} + ${anion.count} ${anion.ion.label}`;
}
function findIonByLabel(label, pool) {
    const wanted = label.trim().toLowerCase().replace(/\s+/g, '');
    for (const ion of pool) {
        if (ion.label.toLowerCase() === wanted)
            return ion;
        // Also accept the plain chemical spelling ("PO43-", "PO4^3-", "SO42-").
        const compact = ion.label.toLowerCase().replace('^', '');
        if (compact === wanted)
            return ion;
    }
    return null;
}
/** Parse an explicit ion override into a candidate. */
function candidateFromLabel(label, pool, side) {
    const ion = findIonByLabel(label, pool);
    if (ion === null) {
        return `${side} "${label}" is not in the recognised ion table (see the README for the full list)`;
    }
    const composition = compositionOf(ion.formula);
    if (composition === null)
        return `${side} "${label}" has an unparseable composition`;
    return { ion, count: 1, composition };
}
function pow(base, exponent) {
    return Math.pow(base, exponent);
}
function formatNumber(value) {
    if (!Number.isFinite(value))
        return String(value);
    if (value === 0)
        return '0';
    const abs = Math.abs(value);
    if (abs >= 1e-3 && abs < 1e5) {
        const fixed = roundOut(value, 6);
        return String(fixed);
    }
    // Compact scientific form for very small/large constants (Ksp ranges widely).
    const exponent = Math.floor(Math.log10(abs));
    const mantissa = roundOut(value / pow(10, exponent), 6);
    return `${mantissa}×10^${exponent}`;
}
/**
 * Solve `(a·s + present)^a · (b·s)^b = ksp` for s ≥ 0 by bisection.
 * `present` is the pre-existing concentration of one of the two ions.
 */
function solveWithCommonIon(ksp, a, b, cationIsCommon, present) {
    const f = (s) => {
        const cat = cationIsCommon ? a * s + present : a * s;
        const an = cationIsCommon ? b * s : b * s + present;
        return pow(cat, a) * pow(an, b) - ksp;
    };
    let lo = 0;
    // Upper bound: the plain-water answer, or the concentration that makes the
    // common-ion side alone exceed Ksp; doubling until f(hi) > 0 is enough.
    let hi = Math.max(1e-12, pow(ksp, 1 / (a + b)));
    for (let i = 0; i < 2000 && f(hi) < 0; i += 1)
        hi *= 2;
    if (!(f(hi) >= 0))
        return null;
    for (let i = 0; i < 400; i += 1) {
        const mid = (lo + hi) / 2;
        if (f(mid) < 0)
            lo = mid;
        else
            hi = mid;
    }
    const s = (lo + hi) / 2;
    const cat = cationIsCommon ? a * s + present : a * s;
    const an = cationIsCommon ? b * s : b * s + present;
    const residual = (pow(cat, a) * pow(an, b) - ksp) / ksp;
    return { s, residual };
}
/**
 * Ksp → molar solubility. See the module header for the grouping rules.
 */
export function kspSolubility(input) {
    const raw = typeof input.formula === 'string' ? input.formula.trim() : '';
    if (raw.length === 0)
        return { ok: false, error: 'formula is empty' };
    const ksp = input.ksp;
    if (typeof ksp !== 'number' || !Number.isFinite(ksp))
        return { ok: false, error: 'ksp must be a finite number' };
    if (ksp <= 0)
        return { ok: false, error: `ksp must be > 0 (got ${ksp}); a non-positive value is not a solubility product` };
    const parse = parseFormula(raw);
    if (!parse.ok)
        return { ok: false, error: parse.error };
    if (parse.hydrate !== null) {
        return {
            ok: false,
            error: `"${raw}" is a hydrate: a Ksp describes the anhydrous solid, so drop the hydrate part (e.g. use CuSO4, not CuSO4·5H2O)`,
        };
    }
    if (parse.coefficient !== null) {
        return { ok: false, error: `"${raw}" starts with a coefficient (${parse.coefficient}); give one formula unit` };
    }
    if (parse.charge !== null) {
        return { ok: false, error: `"${raw}" carries a charge (${parse.charge}); give the neutral salt formula` };
    }
    if (parse.groups.length < 2) {
        return { ok: false, error: `"${raw}" has a single group; a salt needs a cation and an anion` };
    }
    // Enumerate candidate splits: cation half = groups[0..k), anion half = groups[k..).
    const candidates = [];
    for (let k = 1; k < parse.groups.length; k += 1) {
        const cationHalf = halfTotals(parse, 0, k);
        const anionHalf = halfTotals(parse, k, parse.groups.length);
        for (const entry of ION_TABLE) {
            if (entry.ion.charge <= 0)
                continue;
            const a = multipleOf(cationHalf, entry.composition);
            if (a === null)
                continue;
            for (const other of ION_TABLE) {
                if (other.ion.charge >= 0)
                    continue;
                const b = multipleOf(anionHalf, other.composition);
                if (b === null)
                    continue;
                if (a * entry.ion.charge + b * other.ion.charge !== 0)
                    continue;
                candidates.push({
                    cation: { ion: entry.ion, count: a, composition: entry.composition },
                    anion: { ion: other.ion, count: b, composition: other.composition },
                });
            }
        }
    }
    let chosen = null;
    const forcedCation = input.cation !== undefined && input.cation !== '' ? input.cation : null;
    const forcedAnion = input.anion !== undefined && input.anion !== '' ? input.anion : null;
    if (forcedCation !== null || forcedAnion !== null) {
        // An override must still reproduce the formula: the composition and the
        // stoichiometry both have to match, otherwise the caller asked for a salt
        // that is not the one in `formula`.
        const cationHalf = halfTotals(parse, 0, parse.groups.length);
        const allElements = cationHalf;
        const catCandidate = forcedCation !== null ? candidateFromLabel(forcedCation, CATIONS, 'cation') : null;
        if (typeof catCandidate === 'string')
            return { ok: false, error: catCandidate };
        const anCandidate = forcedAnion !== null ? candidateFromLabel(forcedAnion, ANIONS, 'anion') : null;
        if (typeof anCandidate === 'string')
            return { ok: false, error: anCandidate };
        const match = candidates.filter((c) => {
            if (catCandidate !== null && c.cation.ion.label !== catCandidate.ion.label)
                return false;
            if (anCandidate !== null && c.anion.ion.label !== anCandidate.ion.label)
                return false;
            return true;
        });
        if (match.length === 0) {
            const ions = allElements.size;
            return {
                ok: false,
                error: `no ion grouping of "${raw}" matches the requested ions (${ions} element${ions === 1 ? '' : 's'} in the formula); check the charges`,
            };
        }
        if (match.length > 1) {
            return {
                ok: false,
                ambiguous: true,
                error: `several groupings of "${raw}" match the requested ions; give both cation and anion`,
                candidates: match.map((c) => describe(c.cation, c.anion)),
            };
        }
        chosen = match[0];
    }
    else if (candidates.length === 1) {
        chosen = candidates[0];
    }
    else if (candidates.length === 0) {
        return {
            ok: false,
            error: `"${raw}" does not split into two recognised ions; pass cation and anion explicitly (e.g. cation "Fe^3+", anion "OH^-")`,
        };
    }
    else {
        return {
            ok: false,
            ambiguous: true,
            error: `"${raw}" matches ${candidates.length} ion groupings and they do not give the same answer; pass cation and anion to choose`,
            candidates: candidates.map((c) => describe(c.cation, c.anion)),
        };
    }
    const a = chosen.cation.count;
    const b = chosen.anion.count;
    const factor = pow(a, a) * pow(b, b);
    const order = a + b;
    const molarSolubility = pow(ksp / factor, 1 / order);
    const mass = formulaMass(parse);
    const molarMass = mass.molarMass;
    const massSolubility = molarSolubility * molarMass;
    const notes = [];
    if (molarSolubility > 0.1) {
        notes.push(`molar solubility is ${formatNumber(molarSolubility)} mol/L: at this concentration activity coefficients matter, `
            + 'so the ideal-solution Ksp treatment is only roughly right');
    }
    notes.push('Ksp is used exactly as supplied; it is temperature dependent and this tool does not interpolate it');
    let commonIon;
    if (input.commonIon !== undefined && input.commonIon !== '') {
        const label = input.commonIon.trim();
        const isCation = findIonByLabel(label, CATIONS) !== null;
        const isAnion = findIonByLabel(label, ANIONS) !== null;
        if (!isCation && !isAnion) {
            return { ok: false, error: `common ion "${label}" is not in the recognised ion table` };
        }
        const wanted = isCation ? chosen.cation.ion.label : chosen.anion.ion.label;
        const matches = findIonByLabel(label, isCation ? CATIONS : ANIONS);
        if (matches === null || matches.label !== wanted) {
            return {
                ok: false,
                error: `common ion "${label}" is not one of the ions of ${raw} (this salt gives ${chosen.cation.ion.label} and ${chosen.anion.ion.label})`,
            };
        }
        const present = input.commonIonConcentration;
        if (typeof present !== 'number' || !Number.isFinite(present)) {
            return { ok: false, error: 'commonIonConcentration must be a finite number (mol/L) when commonIon is given' };
        }
        if (present <= 0)
            return { ok: false, error: `commonIonConcentration must be > 0 (got ${present})` };
        const solved = solveWithCommonIon(ksp, a, b, isCation, present);
        if (solved === null) {
            return { ok: false, error: 'common-ion solve did not converge; no answer is better than a wrong one' };
        }
        if (!(Math.abs(solved.residual) <= 1e-9)) {
            return {
                ok: false,
                error: `common-ion solution failed its own back-substitution (relative residual ${solved.residual}); refusing to report it`,
            };
        }
        commonIon = {
            symbol: wanted,
            concentration: present,
            molarSolubility: solved.s,
            massSolubility: solved.s * molarMass,
            reductionFactor: solved.s > 0 ? molarSolubility / solved.s : Number.POSITIVE_INFINITY,
            residual: solved.residual,
        };
    }
    const catExpr = a === 1 ? `[${chosen.cation.ion.label}]` : `[${chosen.cation.ion.label}]^${a}`;
    const anExpr = b === 1 ? `[${chosen.anion.ion.label}]` : `[${chosen.anion.ion.label}]^${b}`;
    const catS = a === 1 ? 's' : `(${a}s)`;
    const anS = b === 1 ? 's' : `(${b}s)`;
    const sPart = order === 1 ? 's' : `s^${order}`;
    const factorPart = factor === 1 ? '' : `${factor}·`;
    const expression = `Ksp = ${catExpr}·${anExpr} = ${catS}·${anS} = ${factorPart}${sPart}`;
    const result = {
        ok: true,
        formula: parse.normalized,
        ksp,
        pKsp: roundOut(-Math.log10(ksp), 10),
        cation: {
            symbol: chosen.cation.ion.label,
            charge: chosen.cation.ion.charge,
            count: a,
            concentration: a * molarSolubility,
        },
        anion: {
            symbol: chosen.anion.ion.label,
            charge: chosen.anion.ion.charge,
            count: b,
            concentration: b * molarSolubility,
        },
        stoichiometricFactor: factor,
        ionCount: order,
        expression,
        molarMass: roundOut(molarMass, 6),
        molarSolubility,
        massSolubility,
        chargeBalance: a * chosen.cation.ion.charge + b * chosen.anion.ion.charge,
        notes,
    };
    if (commonIon !== undefined) {
        result.commonIon = commonIon;
    }
    return result;
}
//# sourceMappingURL=ksp.js.map