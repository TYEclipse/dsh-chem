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
import { parseFormula } from "./core.js";
const F_ZERO = { n: 0n, d: 1n };
const F_ONE = { n: 1n, d: 1n };
function bgcd(a, b) {
    let x = a < 0n ? -a : a;
    let y = b < 0n ? -b : b;
    while (y !== 0n) {
        const t = x % y;
        x = y;
        y = t;
    }
    return x;
}
function frac(n, d) {
    if (d === 0n)
        throw new Error('internal: rational with zero denominator');
    let num = n;
    let den = d;
    if (den < 0n) {
        num = -num;
        den = -den;
    }
    const g = bgcd(num, den);
    return g === 0n ? { n: 0n, d: 1n } : { n: num / g, d: den / g };
}
const fAdd = (a, b) => frac(a.n * b.d + b.n * a.d, a.d * b.d);
const fSub = (a, b) => frac(a.n * b.d - b.n * a.d, a.d * b.d);
const fMul = (a, b) => frac(a.n * b.n, a.d * b.d);
const fDiv = (a, b) => frac(a.n * b.d, a.d * b.n);
const fIsZero = (a) => a.n === 0n;
/* ------------------------------------------------------------------ */
/* input parsing                                                       */
/* ------------------------------------------------------------------ */
const ARROW_RE = /-{1,2}>|→|⇌|=>|==>|=/g;
/** A `+` is a species separator under the documented rule (see module docs). */
function isSeparator(text, index) {
    const prev = index === 0 ? undefined : text[index - 1];
    const next = index + 1 >= text.length ? undefined : text[index + 1];
    if (prev === undefined || next === undefined)
        return false;
    if (/\s/.test(prev) || /\s/.test(next)) {
        // Spaced form: "+" with whitespace on both sides only.
        return /\s/.test(prev) && /\s/.test(next);
    }
    return /[A-Za-z)0-9]/.test(prev) && /[A-Z(]/.test(next);
}
/** Split one side of the arrow into species tokens (leading coefficients stripped). */
function splitSide(side) {
    const text = side.trim();
    if (text.length === 0)
        return { error: 'one side of the arrow is empty' };
    const parts = [];
    let current = '';
    for (let i = 0; i < text.length; i += 1) {
        const ch = text[i];
        if (ch === '+' && isSeparator(text, i)) {
            parts.push(current);
            current = '';
        }
        else {
            current += ch;
        }
    }
    parts.push(current);
    const cleaned = [];
    for (const raw of parts) {
        const token = raw.trim().replace(/^\d+/, '');
        if (token.length === 0)
            return { error: `empty species near "${raw.trim()}"` };
        cleaned.push(token);
    }
    return cleaned;
}
/** Parse one species token into counts + charge. */
function parseSpecies(token) {
    const parsed = parseFormula(token);
    if (!parsed.ok)
        return { error: parsed.error };
    const counts = new Map();
    for (const entry of parsed.totals)
        counts.set(entry.symbol, (counts.get(entry.symbol) ?? 0) + entry.count);
    if (parsed.hydrate !== null) {
        for (const entry of parsed.hydrate)
            counts.set(entry.symbol, (counts.get(entry.symbol) ?? 0) + entry.count);
    }
    let charge = 0;
    if (parsed.charge !== null) {
        const sign = parsed.charge.endsWith('-') ? -1 : 1;
        const magnitude = parsed.charge.slice(0, -1);
        charge = sign * (magnitude.length === 0 ? 1 : Number(magnitude));
    }
    return { formula: parsed.normalized, counts, charge };
}
function toRref(matrix, ncols) {
    const rows = matrix.map((row) => row.slice());
    const pivots = [];
    let r = 0;
    for (let c = 0; c < ncols && r < rows.length; c += 1) {
        let pick = -1;
        for (let i = r; i < rows.length; i += 1) {
            if (!fIsZero(rows[i][c])) {
                pick = i;
                break;
            }
        }
        if (pick === -1)
            continue;
        const tmp = rows[r];
        rows[r] = rows[pick];
        rows[pick] = tmp;
        const pivot = rows[r][c];
        for (let j = 0; j < ncols; j += 1)
            rows[r][j] = fDiv(rows[r][j], pivot);
        for (let i = 0; i < rows.length; i += 1) {
            if (i === r)
                continue;
            const factor = rows[i][c];
            if (fIsZero(factor))
                continue;
            for (let j = 0; j < ncols; j += 1) {
                rows[i][j] = fSub(rows[i][j], fMul(factor, rows[r][j]));
            }
        }
        pivots.push(c);
        r += 1;
    }
    return { rows: rows.slice(0, pivots.length), pivots };
}
/** Basis of the nullspace (kernel) of the coefficient matrix. */
function nullspace(matrix, ncols) {
    const { rows, pivots } = toRref(matrix, ncols);
    const isPivot = new Set(pivots);
    const basis = [];
    for (let free = 0; free < ncols; free += 1) {
        if (isPivot.has(free))
            continue;
        const vector = Array.from({ length: ncols }, () => F_ZERO);
        vector[free] = F_ONE;
        rows.forEach((row, index) => {
            const pivotCol = pivots[index];
            vector[pivotCol] = fSub(F_ZERO, row[free]);
        });
        basis.push(vector);
    }
    return basis;
}
/** Clear denominators of a rational vector → integer vector (reduced, sign-normalized). */
function toIntegerVector(vector) {
    let lcm = 1n;
    for (const value of vector) {
        lcm = (lcm / bgcd(lcm, value.d)) * value.d;
    }
    const ints = vector.map((value) => value.n * (lcm / value.d));
    return normalizeIntegers(ints);
}
function normalizeIntegers(ints) {
    let g = 0n;
    for (const value of ints)
        g = bgcd(g, value);
    let out = g === 0n ? ints : ints.map((value) => value / g);
    const firstNonZero = out.find((value) => value !== 0n);
    if (firstNonZero !== undefined && firstNonZero < 0n)
        out = out.map((value) => -value);
    return out;
}
const MAX_COEFFICIENT = 1000000n;
function asPositiveIntegers(ints) {
    const out = [];
    for (const value of ints) {
        if (value <= 0n || value > MAX_COEFFICIENT)
            return null;
        out.push(Number(value));
    }
    return out;
}
/** Smallest-sum integer combination of the basis that is strictly positive. */
function pickSmallestPositive(basis, ncols) {
    const radius = 6;
    const seen = new Set();
    let best = null;
    let bestScore = Number.POSITIVE_INFINITY;
    const combos = (radius * 2 + 1) ** basis.length;
    if (combos > 250_000)
        return null;
    const walk = (index, acc) => {
        if (index === basis.length) {
            const ints = toIntegerVector(acc);
            const key = ints.join(',');
            if (seen.has(key) || ints.every((value) => value === 0n))
                return;
            seen.add(key);
            const positive = asPositiveIntegers(ints);
            if (positive === null)
                return;
            const score = positive.reduce((sum, value) => sum + value, 0);
            if (score < bestScore) {
                bestScore = score;
                best = positive;
            }
            return;
        }
        for (let k = -radius; k <= radius; k += 1) {
            const scaled = basis[index].map((value) => fMul(frac(BigInt(k), 1n), value));
            walk(index + 1, acc.map((value, pos) => fAdd(value, scaled[pos])));
        }
    };
    walk(0, Array.from({ length: ncols }, () => F_ZERO));
    return best;
}
/* ------------------------------------------------------------------ */
/* public entry point                                                  */
/* ------------------------------------------------------------------ */
/** Balance `equation`, e.g. `"C3H8 + O2 -> CO2 + H2O"`. */
export function balanceEquation(equation) {
    const raw = equation ?? '';
    if (typeof raw !== 'string' || raw.trim().length === 0) {
        return { ok: false, kind: 'parse', error: 'equation is empty' };
    }
    if (raw.length > 400) {
        return { ok: false, kind: 'parse', error: `equation is too long (${raw.length} chars; max 400)` };
    }
    const arrows = raw.match(ARROW_RE) ?? [];
    if (arrows.length !== 1) {
        return {
            ok: false,
            kind: 'parse',
            error: `expected exactly one reaction arrow (-> or =), found ${arrows.length}`,
        };
    }
    const arrowIndex = raw.search(ARROW_RE);
    const arrowText = arrows[0];
    const leftRaw = raw.slice(0, arrowIndex);
    const rightRaw = raw.slice(arrowIndex + arrowText.length);
    const leftTokens = splitSide(leftRaw);
    if (!Array.isArray(leftTokens))
        return { ok: false, kind: 'parse', error: `left side: ${leftTokens.error}` };
    const rightTokens = splitSide(rightRaw);
    if (!Array.isArray(rightTokens))
        return { ok: false, kind: 'parse', error: `right side: ${rightTokens.error}` };
    const reactants = [];
    for (const token of leftTokens) {
        const species = parseSpecies(token);
        if ('error' in species)
            return { ok: false, kind: 'parse', error: `reactant "${token}": ${species.error}` };
        reactants.push(species);
    }
    const products = [];
    for (const token of rightTokens) {
        const species = parseSpecies(token);
        if ('error' in species)
            return { ok: false, kind: 'parse', error: `product "${token}": ${species.error}` };
        products.push(species);
    }
    const ncols = reactants.length + products.length;
    if (ncols < 2) {
        return { ok: false, kind: 'parse', error: 'at least two species are required (one reactant and one product)' };
    }
    const elements = [];
    const seenElements = new Set();
    for (const species of [...reactants, ...products]) {
        for (const symbol of species.counts.keys()) {
            if (!seenElements.has(symbol)) {
                seenElements.add(symbol);
                elements.push(symbol);
            }
        }
    }
    if (elements.length === 0) {
        return { ok: false, kind: 'parse', error: 'no elements found on either side' };
    }
    const anyCharge = [...reactants, ...products].some((species) => species.charge !== 0);
    const matrix = elements.map((symbol) => [...reactants, ...products].map((species, index) => {
        const count = BigInt(species.counts.get(symbol) ?? 0);
        return frac(index < reactants.length ? count : -count, 1n);
    }));
    if (anyCharge) {
        matrix.push([...reactants, ...products].map((species, index) => {
            const charge = BigInt(species.charge);
            return frac(index < reactants.length ? charge : -charge, 1n);
        }));
    }
    const basis = nullspace(matrix, ncols);
    if (basis.length === 0) {
        return {
            ok: false,
            kind: 'inconsistent',
            error: 'no nonzero solution: the elements listed cannot be conserved with these species',
        };
    }
    let solution;
    let unique = true;
    if (basis.length === 1) {
        solution = asPositiveIntegers(toIntegerVector(basis[0]));
    }
    else {
        unique = false;
        solution = pickSmallestPositive(basis, ncols);
    }
    if (solution === null) {
        if (basis.length > 1) {
            return {
                ok: false,
                kind: 'too-ambitious',
                error: `more than one independent balance exists (${basis.length} free parameters): `
                    + 'all species cannot take part in a single equation — split it into separate reactions',
            };
        }
        return {
            ok: false,
            kind: 'no-positive',
            error: 'no strictly positive integer solution: at least one species would need coefficient 0 '
                + '(it does not take part in this reaction) or the equation cannot be balanced',
        };
    }
    // Defensive self-check: recount atoms from the returned coefficients.
    const reported = [];
    for (const symbol of elements) {
        let left = 0;
        let right = 0;
        reactants.forEach((species, index) => {
            left += solution[index] * (species.counts.get(symbol) ?? 0);
        });
        products.forEach((species, index) => {
            right += solution[reactants.length + index] * (species.counts.get(symbol) ?? 0);
        });
        if (left !== right) {
            return { ok: false, kind: 'inconsistent', error: `internal check failed: element ${symbol} is not conserved` };
        }
        reported.push({ symbol, left, right });
    }
    return {
        ok: true,
        reactants: reactants.map((species, index) => ({ formula: species.formula, coefficient: solution[index] })),
        products: products.map((species, index) => ({
            formula: species.formula,
            coefficient: solution[reactants.length + index],
        })),
        elements: reported,
        unique,
        chargeBalanced: anyCharge,
    };
}
/** Render a balanced result as `2 Fe + 3 O2 -> Fe2O3`-style text. */
export function renderBalancedEquation(reactants, products) {
    const side = (terms) => terms.map((term) => `${term.coefficient > 1 ? `${term.coefficient} ` : ''}${term.formula}`).join(' + ');
    return `${side(reactants)} -> ${side(products)}`;
}
//# sourceMappingURL=balance.js.map