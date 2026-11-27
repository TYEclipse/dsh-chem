#!/usr/bin/env python3
"""
Independent anchor oracle for dsh-chem tests.

Deliberately independent implementation (own periodic-table dict typed from the
IUPAC standard-atomic-weight table, own tokenizer) used ONLY to generate
expected values for the vitest suite. Not shipped in the plugin.

Charge disambiguation rule (mirrors the TypeScript parser, documented in the
README): digits immediately before a trailing +/- belong to the IONIC CHARGE
only when every element parsed so far is the same single element
("Fe3+" -> Fe with 3+, "S2-" -> S with 2-, "Ca2+" -> Ca with 2+). Otherwise
they are a subscript and the bare sign is a 1-unit charge ("NH4+" -> N + H4
with +, "NO3-" -> N + O3 with -). The explicit form "SO4^2-" is always the
charge. "(O2)+" escapes the rule for subscript-before-charge.

Cross-validation: anchors must match publicly published textbook values
(H2O 18.015, H2SO4 98.08, NaCl 58.44, Ca(OH)2 74.09, CuSO4*5H2O 249.68,
C6H12O6 180.16, KMnO4 158.03, NH3 17.03).
"""
import itertools
import math, re, sys
from fractions import Fraction

ATOMIC = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.95,
    'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41,
    'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29,
    'Cs': 132.91, 'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
    'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.50,
    'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05, 'Lu': 174.97, 'Hf': 178.49,
    'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08,
    'Au': 196.97, 'Hg': 200.59, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98,
    'Po': 209.0, 'At': 210.0, 'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0,
    'Th': 232.04, 'Pa': 231.04, 'U': 238.03, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0,
    'Cm': 247.0, 'Bk': 247.0, 'Cf': 251.0, 'Es': 252.0, 'Fm': 257.0, 'Md': 258.0,
    'No': 259.0, 'Lr': 266.0, 'Rf': 267.0, 'Db': 268.0, 'Sg': 269.0, 'Bh': 270.0,
    'Hs': 277.0, 'Mt': 278.0, 'Ds': 281.0, 'Rg': 282.0, 'Cn': 285.0, 'Nh': 286.0,
    'Fl': 289.0, 'Mc': 290.0, 'Lv': 293.0, 'Ts': 294.0, 'Og': 294.0,
}
AVOGADRO = 6.02214076e23
SINGLE_LOWERCASE_OK = set('bcfhiknopsuvwy')  # lowercase letters that ARE symbols

def parse_segment(seg, original):
    tokens = re.findall(r'[A-Z][a-z]?|[a-z]|\d+|\(|\)|[+-]', seg)
    counts = {}
    stack = [counts]
    charge = None
    i, n = 0, len(tokens)
    while i < n:
        t = tokens[i]
        if t == '(':
            stack.append({}); i += 1
            continue
        if t == ')':
            inner = stack.pop(); i += 1
            if len(stack) == 0:
                raise ValueError(f'unbalanced parentheses in {original}')
            top = stack[-1]
            if i < n and tokens[i].isdigit() and i + 1 < n and tokens[i + 1] in '+-' and len(set(inner.keys())) == 1:
                # "(Fe)3+": digits+sign after a single-element group = charge.
                for k, v in inner.items():
                    top[k] = top.get(k, 0) + v
                charge = tokens[i] + tokens[i + 1]; i += 2
                if i != n: raise ValueError(f'content after charge in {original}')
                break
            mult = 1
            if i < n and tokens[i].isdigit():
                mult = int(tokens[i]); i += 1
                if mult == 0:
                    raise ValueError(f'zero multiplicity in {original}')
            for k, v in inner.items():
                top[k] = top.get(k, 0) + v * mult
            if i < n and tokens[i] in '+-':
                charge = tokens[i]; i += 1
                if i != n: raise ValueError(f'content after charge in {original}')
                break
            continue
        if t.isdigit():
            raise ValueError(f'stray number {t} in {original}')
        sym = t
        if len(sym) == 1 and sym.islower():
            if sym in SINGLE_LOWERCASE_OK:
                sym = sym.upper()
            else:
                raise ValueError(f'unknown element symbol "{sym}" in {original}')
        if sym not in ATOMIC:
            raise ValueError(f'unknown element {sym} in {original}')
        i += 1
        mult = 1
        if i < n and tokens[i].isdigit():
            digits = tokens[i]
            distinct = set()
            for layer in stack:
                for k in layer:
                    distinct.add(k)
            distinct.add(sym)
            if i + 1 < n and tokens[i + 1] in '+-' and len(distinct) == 1 and len(stack) == 1:
                # "Fe3+": digits+sign after the only element = charge magnitude.
                top = stack[-1]
                top[sym] = top.get(sym, 0) + 1
                charge = digits + tokens[i + 1]; i += 2
                if i != n: raise ValueError(f'content after charge in {original}')
                break
            mult = int(digits); i += 1
            if mult == 0:
                raise ValueError(f'zero multiplicity in {original}')
        top = stack[-1]
        top[sym] = top.get(sym, 0) + mult
        if i < n and tokens[i] in '+-':
            # "Na+" / "NH4+": bare trailing sign = 1-unit charge.
            charge = tokens[i]; i += 1
            if i != n: raise ValueError(f'content after charge in {original}')
            break
    if len(stack) != 1:
        raise ValueError(f'unbalanced parentheses in {original}')
    return counts, charge

def parse_formula(formula):
    s = formula.strip().replace('·', '.').replace('*', '.').replace('•', '.')
    if s == '':
        raise ValueError('formula is empty')
    # Explicit ^ charge suffix (unambiguous): SO4^2-, Fe^3+
    m = re.search(r'\^(\d*)([+-])\s*$', s)
    explicit_charge = None
    if m:
        explicit_charge = (m.group(1) if m.group(1) else '') + m.group(2)
        s = s[:m.start()].rstrip('.')
    parts = s.split('.')
    if len(parts) > 2:
        raise ValueError(f'too many hydrate separators in {formula}')
    main_seg = parts[0]
    coeff = None
    m = re.match(r'^(\d+)', main_seg)
    if m and len(m.group(1)) < len(main_seg):
        coeff = int(m.group(1))
        main_seg = main_seg[len(m.group(1)):]
    counts, charge = parse_segment(main_seg, formula)
    if explicit_charge is not None:
        if charge is not None:
            raise ValueError(f'double charge in {formula}')
        charge = explicit_charge
    hydrate = None
    if len(parts) == 2:
        hseg = parts[1]
        hydrate_mult = 1
        m = re.match(r'^(\d+)', hseg)
        if m:
            hydrate_mult = int(m.group(1))
            hseg = hseg[len(m.group(1)):]
        hcounts, hcharge = parse_segment(hseg, formula)
        if hcharge:
            if charge is not None:
                raise ValueError(f'double charge in {formula}')
            charge = hcharge
        hydrate = {k: v * hydrate_mult for k, v in hcounts.items()}
    total = dict(counts)
    if hydrate:
        for k, v in hydrate.items():
            total[k] = total.get(k, 0) + v
    mass = sum(ATOMIC[k] * v for k, v in total.items())
    return mass, total, hydrate, charge, coeff

TEST_FORMULAS = [
    'H2O', 'H2SO4', 'NaCl', 'Ca(OH)2', 'CuSO4·5H2O', 'C6H12O6', 'CH3CH2OH',
    'Fe2(SO4)3', 'KMnO4', 'NH3', 'KAl(SO4)2·12H2O', 'NaHCO3', 'CO2', 'CoCl2',
    'Fe3+', 'Fe2+', 'SO4^2-', 'NH4+', 'NO3-', 'PO4^3-', '2H2O', 'CaCO3',
    'AgNO3', 'C2H5OH', 'CH3COOH', 'Mg(OH)2', '(NH4)2SO4', 'U', 'Tc', 'h2o',
    'o2', 'CaCl2·2H2O', 'S2-', 'O2-', 'N3-', 'MnO4-', 'O2^+', 'CuSO4',
    '(CH3)2CO', 'CH3COCH3', 'Na+', 'OH-', 'Cl-', 'K+', 'Mg2+', 'Al3+',
]
print('# --- molar mass anchors (independent parser) ---')
for f in TEST_FORMULAS:
    try:
        mass, total, hydrate, charge, coeff = parse_formula(f)
        items = ','.join(f'{k}:{v}' for k, v in sorted(total.items()))
        print(f'{f}\t{mass:.9f}\t{items}\thydrate={hydrate}\tcharge={charge}\tcoeff={coeff}')
    except ValueError as e:
        print(f'{f}\tERROR\t{e}')

print('# --- textbook cross-validation (published values) ---')
TEXTBOOK = {'H2O': 18.015, 'H2SO4': 98.08, 'NaCl': 58.44, 'Ca(OH)2': 74.09,
            'CuSO4·5H2O': 249.68, 'C6H12O6': 180.16, 'KMnO4': 158.03, 'NH3': 17.03}
all_ok = True
for f, expected in TEXTBOOK.items():
    mass, _, _, _, _ = parse_formula(f)
    delta = abs(mass - expected)
    ok = delta <= 0.011
    all_ok = all_ok and ok
    print(f'{f}: computed={mass:.3f} textbook={expected} delta={delta:.4f} {"OK" if ok else "FAIL"}')
print('TEXTBOOK-ALL-OK' if all_ok else 'TEXTBOOK-HAS-FAILURES')

print('# --- conversion anchors ---')
M_H2SO4 = parse_formula('H2SO4')[0]
print(f'98g_H2SO4_to_mol\t{98 / M_H2SO4:.12f}')
print(f'0.5mol_H2SO4_to_g\t{0.5 * M_H2SO4:.12f}')
print(f'10g_NaCl_to_particles\t{(10 / parse_formula("NaCl")[0]) * AVOGADRO:.6e}')
print(f'1e24_particles_water_to_g\t{(1e24 / AVOGADRO) * parse_formula("H2O")[0]:.9f}')
print(f'2.5mol_to_particles\t{2.5 * AVOGADRO:.6e}')
print(f'1.5e23_particles_to_mol\t{1.5e23 / AVOGADRO:.12f}')

print('# --- solution anchors ---')
print(f'molarity_NaCl_58.44g_250ml\t{(58.44 / parse_formula("NaCl")[0]) / 0.25:.12f}')
print(f'prepare_0.1M_NaOH_500ml\t{0.1 * 0.5 * parse_formula("NaOH")[0]:.9f}')
print(f'dilute_v1\t{(0.5 * 250) / 2.0:.9f}')
print(f'dilute_c2\t{(2.0 * 100) / 500:.9f}')
print(f'dilute_v2\t{(2.0 * 100) / 0.5:.9f}')
print(f'dilute_c1\t{(0.5 * 250) / 10:.9f}')
print(f'pct_wv\t{(5.0 / 200.0) * 100:.9f}')
print(f'pct_vv\t{(30.0 / 1000.0) * 100:.9f}')
print(f'pct_mass\t{(2.5 / 100.0) * 100:.9f}')

print('# --- pH anchors ---')
print(f'ph_from_h_1e-7\t{-math.log10(1e-7):.12f}')
print(f'ph_from_h_0.01\t{-math.log10(0.01):.12f}')
print(f'ph_from_oh_1e-4\t{14 - (-math.log10(1e-4)):.12f}')
print(f'ph_from_ph_3_h\t{math.pow(10, -3):.6e}')
print(f'ph_from_poh_2_oh\t{math.pow(10, -2):.6e}')
print(f'buffer_4.76_0.1_0.2\t{4.76 + math.log10(0.1 / 0.2):.12f}')
print(f'buffer_9.25_0.05_0.05\t{9.25 + math.log10(0.05 / 0.05):.12f}')
print(f'buffer_textbook_check\t{4.76 + math.log10(0.1 / 0.2):.4f} (expect ~4.46)')

# ---------------------------------------------------------------------------
# balance_equation anchors (independent Fraction-based nullspace solver)
#
# Independent implementation: Python's fractions.Fraction + own species
# tokenizer + own RREF. Deterministic selection rule for the multi-solution
# case: among all strictly positive integer solutions reachable from the
# nullspace basis with small integer combinations, take the SMALLEST TOTAL SUM
# (ties broken by enumeration order = lexicographic in the basis coefficients).
print()
print('# --- balance_equation anchors (Fraction nullspace) ---')

def _is_sep(text, i):
    prev = text[i - 1] if i > 0 else None
    nxt = text[i + 1] if i + 1 < len(text) else None
    if prev is None or nxt is None:
        return False
    if prev.isspace() or nxt.isspace():
        return prev.isspace() and nxt.isspace()
    return (prev.isalnum() or prev == ')') and (nxt.isupper() or nxt == '(')


def _split_species(side):
    text = side.strip()
    if not text:
        raise ValueError('one side of the arrow is empty')
    parts, cur = [], ''
    for i, ch in enumerate(text):
        if ch == '+' and _is_sep(text, i):
            parts.append(cur)
            cur = ''
        else:
            cur += ch
    parts.append(cur)
    out = []
    for raw in parts:
        token = re.sub(r'^\d+', '', raw.strip())
        if not token:
            raise ValueError(f'empty species near {raw!r}')
        out.append(token)
    return out


def _charge_of(charge):
    if not charge:
        return 0
    sign = -1 if charge.endswith('-') else 1
    mag = charge[:-1]
    return sign * (int(mag) if mag else 1)


def _rref(matrix, ncols):
    rows = [[Fraction(x) for x in row] for row in matrix]
    pivots, r = [], 0
    for c in range(ncols):
        if r >= len(rows):
            break
        pick = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pick is None:
            continue
        rows[r], rows[pick] = rows[pick], rows[r]
        pivot = rows[r][c]
        rows[r] = [v / pivot for v in rows[r]]
        for i in range(len(rows)):
            if i == r:
                continue
            factor = rows[i][c]
            if factor == 0:
                continue
            rows[i] = [a - factor * b for a, b in zip(rows[i], rows[r])]
        pivots.append(c)
        r += 1
    return rows[:len(pivots)], pivots


def _nullspace(matrix, ncols):
    rows, pivots = _rref(matrix, ncols)
    basis = []
    for free in [c for c in range(ncols) if c not in pivots]:
        vec = [Fraction(0)] * ncols
        vec[free] = Fraction(1)
        for i, row in enumerate(rows):
            vec[pivots[i]] = -row[free]
        basis.append(vec)
    return basis


def _integers(vec):
    lcm = 1
    for x in vec:
        lcm = lcm * x.denominator // math.gcd(lcm, x.denominator)
    ints = [int(x * lcm) for x in vec]
    g = 0
    for v in ints:
        g = math.gcd(g, abs(v))
    if g:
        ints = [v // g for v in ints]
    first = next((v for v in ints if v != 0), 0)
    if first < 0:
        ints = [-v for v in ints]
    return ints


ARROW_RE_PY = re.compile(r'-{1,2}>|→|⇌|=>|==>|=')
MAX_COEFF = 1_000_000


def balance_equation(eq):
    """Return (kind, coefficients, unique, elements) for a chemical equation."""
    found = ARROW_RE_PY.findall(eq)
    if len(found) != 1:
        raise ValueError(f'expected exactly one arrow, found {len(found)}')
    at = ARROW_RE_PY.search(eq)
    left = _split_species(eq[:at.start()])
    right = _split_species(eq[at.end():])
    species = []
    for token in left:
        _m, total, _h, charge, _c = parse_formula(token)
        species.append((total, _charge_of(charge)))
    order = len(left)
    for token in right:
        _m, total, _h, charge, _c = parse_formula(token)
        species.append((total, _charge_of(charge)))
    ncols = len(species)
    if ncols < 2:
        raise ValueError('need at least two species')
    elements = []
    for total, _c in species:
        for k in total:
            if k not in elements:
                elements.append(k)
    matrix = []
    for el in elements:
        matrix.append([Fraction(total.get(el, 0)) * (1 if i < order else -1)
                       for i, (total, _c) in enumerate(species)])
    if any(c for _t, c in species):
        matrix.append([Fraction(c) * (1 if i < order else -1) for i, (_t, c) in enumerate(species)])
    basis = _nullspace(matrix, ncols)
    if not basis:
        return ('inconsistent', [], True, elements)
    if len(basis) == 1:
        coeffs, unique = _integers(basis[0]), True
    else:
        unique, coeffs, best_sum = False, None, None
        radius = 6
        for combo in itertools.product(range(-radius, radius + 1), repeat=len(basis)):
            if all(k == 0 for k in combo):
                continue
            acc = [Fraction(0)] * ncols
            for k, base in zip(combo, basis):
                acc = [a + k * b for a, b in zip(acc, base)]
            ints = _integers(acc)
            if all(0 < v <= MAX_COEFF for v in ints):
                total_sum = sum(ints)
                if best_sum is None or total_sum < best_sum:
                    best_sum, coeffs = total_sum, ints
    if coeffs is None or not all(v > 0 for v in coeffs):
        return ('no-positive', [], unique, elements)
    # Verify conservation from the returned coefficients (independent recount).
    for index, el in enumerate(elements):
        lhs = sum(coeffs[i] * species[i][0].get(el, 0) for i in range(order))
        rhs = sum(coeffs[i] * species[i][0].get(el, 0) for i in range(order, ncols))
        if lhs != rhs:
            raise AssertionError(f'{eq}: {el} not conserved ({lhs} != {rhs})')
    if any(c for _t, c in species):
        lhs = sum(coeffs[i] * species[i][1] for i in range(order))
        rhs = sum(coeffs[i] * species[i][1] for i in range(order, ncols))
        if lhs != rhs:
            raise AssertionError(f'{eq}: charge not conserved ({lhs} != {rhs})')
    return ('ok', coeffs, unique, elements)


BALANCE_CASES = [
    # textbook anchors (published coefficients)
    'C3H8 + O2 -> CO2 + H2O',            # 1 5 3 4
    'Fe + O2 -> Fe2O3',                  # 4 3 2
    'KClO3 -> KCl + O2',                 # 2 2 3
    'Fe2O3 + CO -> Fe + CO2',            # 1 3 2 3
    'Al + HCl -> AlCl3 + H2',            # 2 6 2 3
    'C6H12O6 + O2 -> CO2 + H2O',         # 1 6 6 6
    'Mg(OH)2 + H3PO4 -> Mg3(PO4)2 + H2O',  # 3 2 1 6
    'Cu + HNO3 -> Cu(NO3)2 + NO2 + H2O',   # 1 4 1 2 2
    'Ca(OH)2 + H3PO4 -> Ca3(PO4)2 + H2O',  # 3 2 1 6
    'NH3 + O2 -> NO + H2O',              # 4 5 4 6
    'CH4 + O2 -> CO2 + H2O',             # 1 2 1 2
    'Fe + HCl -> FeCl2 + H2',            # 1 2 1 1
    'NaHCO3 -> Na2CO3 + H2O + CO2',      # 2 1 1 1
    'KMnO4 + HCl -> KCl + MnCl2 + Cl2 + H2O',   # 2 16 2 2 5 8
    # compact separator form and ignored leading coefficients
    '2H2 + O2 -> 2H2O',
    'Fe+O2->Fe2O3',
    # ionic / redox (charge row active)
    'Fe2+ + Ce4+ -> Fe3+ + Ce3+',
    'MnO4- + Fe2+ + H+ -> Mn2+ + Fe3+ + H2O',   # 1 5 8 1 5 4
    'Na+ + Cl- -> NaCl',
    # hydration water participates in the balance
    'CuSO4·5H2O -> CuSO4 + H2O',          # 1 1 5
    'H2 + O2 -> H2O',                     # 2 1 2
]
for case in BALANCE_CASES:
    try:
        kind, coeffs, unique, elements = balance_equation(case)
        joined = ','.join(str(c) for c in coeffs)
        print(f'balance\t{case}\t{kind}\t{joined}\tunique={unique}\telements={len(elements)}')
    except (ValueError, AssertionError) as exc:
        print(f'balance\t{case}\tERROR\t{exc}')

print('# --- balance: multi-solution and impossible cases (properties only) ---')
for case in ['Fe + O2 -> Fe2O3 + Fe3O4', 'H2 + O2 -> H2O + N2', 'Na + Cl2 -> NaCl + KCl']:
    try:
        kind, coeffs, unique, _elements = balance_equation(case)
        print(f'balance-edge\t{case}\t{kind}\tsum={sum(coeffs)}\tunique={unique}\t{",".join(str(c) for c in coeffs)}')
    except (ValueError, AssertionError) as exc:
        print(f'balance-edge\t{case}\tERROR\t{exc}')

# ---------------------------------------------------------------------------
# gas_law anchors (ideal gas law, own constant derivation)
print()
print('# --- gas_law anchors ---')
GAS_R_J = 8.31446261815324
GAS_ATM_KPA = 101.325
GAS_R = GAS_R_J / GAS_ATM_KPA
print(f'R_atm_l\t{GAS_R!r}')
print(f'R_from_joules\t{(GAS_R_J / GAS_ATM_KPA)!r}')
print(f'R_published_check\t{GAS_R:.5f} (published 0.08206)')
print(f'molar_volume_stp_273.15K_1atm\t{(GAS_R * 273.15 / 1)!r}')
print(f'molar_volume_published_check\t{(GAS_R * 273.15):.3f} (published 22.414)')
print(f'molar_volume_298.15K_1atm\t{(GAS_R * 298.15 / 1)!r}')
print(f'molar_volume_273.15K_2atm\t{(GAS_R * 273.15 / 2)!r}')
# ideal-gas solves: (P,V,n,T) with one missing
print(f'ideal_pressure_n1_T273.15_V22.4\t{((1 * GAS_R) * 273.15) / 22.4!r}')
print(f'ideal_volume_n2_T300_P1.5\t{((2 * GAS_R) * 300) / 1.5!r}')
print(f'ideal_moles_P1_V1_T298.15\t{(1 * 1) / (GAS_R * 298.15)!r}')
print(f'ideal_temperature_P2_V10_n1\t{(2 * 10) / (1 * GAS_R)!r}')
print(f'ideal_pressure_kpa_n1_T273.15_V22.4\t{(((1 * GAS_R) * 273.15) / 22.4) * GAS_ATM_KPA!r}')
print(f'ideal_volume_ml_n2_T300_P1.5\t{(((2 * GAS_R) * 300) / 1.5) * 1000!r}')
print(f'ideal_temperature_c_P2_V10_n1\t{((2 * 10) / (1 * GAS_R)) - 273.15!r}')
# density / molar mass from density
CO2_MASS = parse_formula('CO2')[0]
N2_MASS = parse_formula('N2')[0]
print(f'co2_molar_mass\t{CO2_MASS!r}')
print(f'n2_molar_mass\t{N2_MASS!r}')
print(f'density_CO2_1atm_273.15K\t{(1 * CO2_MASS) / (GAS_R * 273.15)!r}')
print(f'density_N2_1atm_298.15K\t{(1 * N2_MASS) / (GAS_R * 298.15)!r}')
print(f'density_CO2_2atm_300K\t{(2 * CO2_MASS) / (GAS_R * 300)!r}')
print(f'density_CO2_g_per_ml_1atm_273.15K\t{((1 * CO2_MASS) / (GAS_R * 273.15)) / 1000!r}')
print(f'molar_mass_from_density_1.96_1atm_273.15K\t{(1.96 * GAS_R * 273.15) / 1!r}')
print(f'molar_mass_from_density_1.2506_1atm_273.15K\t{(1.2506 * GAS_R * 273.15) / 1!r}')
print(f'molar_mass_from_density_1.15_1.5atm_300K\t{(1.15 * GAS_R * 300) / 1.5!r}')
print(f'density_published_check_CO2\t{((1 * CO2_MASS) / (GAS_R * 273.15)):.4f} (published 1.96 g/L)')

# ---------------------------------------------------------------------------
# stoichiometry anchors (limiting reagent / extent / yields)
# Independent of the TS implementation: coefficients come from this file's own
# Fraction balancer, molar masses from this file's own periodic table.
print()
print('# --- stoichiometry anchors (extent / limiting reagent / yields) ---')


def _fmass(species):
    return parse_formula(species)[0]


def stoich(eq, amounts, product=None, actual_g=None):
    """(extent, reactant rows, product rows, limiting names, yield tuple)."""
    at = ARROW_RE_PY.search(eq)
    left = _split_species(eq[:at.start()])
    right = _split_species(eq[at.end():])
    kind, coeffs, _u, _e = balance_equation(eq)
    assert kind == 'ok', (eq, kind)
    rc = dict(zip(left, coeffs[:len(left)]))
    pc = dict(zip(right, coeffs[len(left):]))
    moles = {}
    for name, spec in amounts.items():
        if 'moles' in spec:
            moles[name] = spec['moles']
        else:
            moles[name] = spec['grams'] / _fmass(name)
    extent = min(moles[n] / rc[n] for n in left)
    limiting = [n for n in left if abs(moles[n] / rc[n] - extent) <= 1e-9 * max(1.0, extent)]
    reactants = []
    for n in left:
        m = _fmass(n)
        consumed = extent * rc[n]
        left_m = max(0.0, moles[n] - consumed)
        reactants.append((n, rc[n], m, moles[n], moles[n] * m, consumed, left_m, left_m * m))
    products = []
    for n in right:
        m = _fmass(n)
        mol = extent * pc[n]
        products.append((n, pc[n], m, mol, mol * m))
    y = None
    if product is not None:
        row = [p for p in products if p[0] == product][0]
        if actual_g is None:
            y = (row[0], row[3], row[4])
        else:
            y = (row[0], row[3], row[4], actual_g, actual_g / row[4] * 100)
    return extent, reactants, products, limiting, y


STOICH_CASES = [
    ('propane_stoich', 'C3H8 + O2 -> CO2 + H2O', {'C3H8': {'moles': 1.0}, 'O2': {'moles': 5.0}}, None, None),
    ('propane_oxygen_limiting', 'C3H8 + O2 -> CO2 + H2O', {'C3H8': {'moles': 0.5}, 'O2': {'moles': 2.0}}, 'CO2', None),
    ('propane_percent_yield', 'C3H8 + O2 -> CO2 + H2O', {'C3H8': {'moles': 0.5}, 'O2': {'moles': 2.0}}, 'CO2', 44.0),
    ('zinc_acid_grams', 'Zn + 2HCl -> ZnCl2 + H2', {'Zn': {'grams': 6.54}, 'HCl': {'moles': 0.2}}, 'H2', None),
    ('water_colimiting', 'H2 + O2 -> H2O', {'H2': {'moles': 2.0}, 'O2': {'moles': 1.0}}, 'H2O', None),
    ('hydrate_decomposition', 'CuSO4·5H2O -> CuSO4 + H2O', {'CuSO4·5H2O': {'moles': 0.5}}, 'H2O', None),
    ('permanganate_ionic', 'MnO4- + Fe2+ + H+ -> Mn2+ + Fe3+ + H2O',
     {'MnO4-': {'moles': 0.02}, 'Fe2+': {'moles': 0.1}, 'H+': {'moles': 0.2}}, 'Mn2+', None),
]
for slug, eq, amounts, product, actual in STOICH_CASES:
    extent, reactants, products, limiting, y = stoich(eq, amounts, product, actual)
    print(f'stoich.{slug}.extent\t{extent!r}')
    print(f'stoich.{slug}.limiting\t{",".join(limiting)}')
    for n, coef, m, avail, avail_g, consumed, left_m, left_g in reactants:
        print(f'stoich.{slug}.r.{n}\tcoef={coef}\tM={m!r}\tavailMol={avail!r}\tavailG={avail_g!r}'
              f'\tconsumed={consumed!r}\tleftMol={left_m!r}\tleftG={left_g!r}')
    for n, coef, m, mol, g in products:
        print(f'stoich.{slug}.p.{n}\tcoef={coef}\tM={m!r}\tmol={mol!r}\tg={g!r}')
    if y is not None:
        if len(y) == 3:
            print(f'stoich.{slug}.yield\t{y[0]}\tmol={y[1]!r}\tg={y[2]!r}')
        else:
            print(f'stoich.{slug}.yield\t{y[0]}\tmol={y[1]!r}\tg={y[2]!r}\tactualG={y[3]!r}\tpercent={y[4]!r}')
    # Same anchors at the precision the tests assert them (6 significant digits):
    # the coverage hint compares literals, so this keeps every asserted decimal
    # traceable to a line of oracle output instead of "close enough to the repr".
    for n, coef, m, avail, avail_g, consumed, left_m, left_g in reactants:
        print(f'stoich.{slug}.r.{n}.round6\tM={m:.6g}\tavailMol={avail:.6g}\tavailG={avail_g:.6g}'
              f'\tconsumed={consumed:.6g}\tleftMol={left_m:.6g}\tleftG={left_g:.6g}')
    for n, coef, m, mol, g in products:
        print(f'stoich.{slug}.p.{n}.round6\tM={m:.6g}\tmol={mol:.6g}\tg={g:.6g}')
    if y is not None and len(y) == 5:
        print(f'stoich.{slug}.percent.round6\t{y[4]:.6g}')
    for n, coef, m, avail, avail_g, consumed, left_m, left_g in reactants:
        print(f'stoich.{slug}.r.{n}.round12\tavailMol={avail:.12g}\tavailG={avail_g:.12g}\tleftMol={left_m:.12g}')
    if y is not None and len(y) == 5:
        print(f'stoich.{slug}.percent.round12\t{y[4]:.12g}')

# ---------------------------------------------------------------------------
# formula_from_composition anchors (percent composition / masses -> empirical)
# Independent of the TS implementation: own table, own Hill ordering, own
# smallest-multiplier search. Published compositions are the ruler (water
# 11.19 % H / 88.81 % O, glucose 40.00 / 6.71 / 53.29 with M = 180.156...).
print()
print('# --- formula_from_composition anchors ---')


def _moles_from(comp, mode):
    if mode == 'percent':
        total = sum(comp.values())
        if abs(total - 100.0) > 1.0:
            raise SystemExit(f'percent sum {total} out of range')
        grams = dict(comp)
    else:
        grams = dict(comp)
    return {el: g / ATOMIC[el] for el, g in grams.items()}


def _hill(atoms):
    if 'C' in atoms:
        order = ['C'] + (['H'] if 'H' in atoms else []) + sorted(el for el in atoms if el not in ('C', 'H'))
    else:
        order = sorted(atoms)
    out = ''
    for el in order:
        n = atoms[el]
        out += el if n == 1 else f'{el}{n}'
    return out


def empirical(comp, mode, molar_mass=None, tol=0.05, kmax=8):
    moles = _moles_from(comp, mode)
    mn = min(moles.values())
    ratios = {el: m / mn for el, m in moles.items()}
    best_k, best_dev = 1, None
    exact = False
    for k in range(1, kmax + 1):
        dev = max(abs(r * k - round(r * k)) for r in ratios.values())
        if best_dev is None or dev < best_dev - 1e-15:
            best_k, best_dev = k, dev
        if dev <= tol:
            best_k, best_dev, exact = k, dev, True
            break
        # keep the smallest k on ties: only replace on a strictly smaller dev
    atoms = {el: max(1, int(round(r * best_k))) for el, r in ratios.items()}
    mass = sum(atoms[el] * ATOMIC[el] for el in atoms)
    pct = {el: atoms[el] * ATOMIC[el] / mass * 100 for el in atoms}
    rec = dict(atoms=atoms, kmult=best_k, dev=best_dev, exact=exact, formula=_hill(atoms),
               mass=mass, pct=pct, moles=moles, ratios=ratios)
    if molar_mass is not None:
        n = molar_mass / mass
        nearest = max(1, int(round(n)))
        rel = abs(n - nearest) / nearest
        rec.update(nmol=n, nearest=nearest, nrel=rel, nexact=rel <= 0.02,
                   implied=nearest * mass, merr=abs(nearest * mass - molar_mass) / molar_mass * 100,
                   molatoms={el: atoms[el] * nearest for el in atoms},
                   molformula=_hill({el: atoms[el] * nearest for el in atoms}))
    return rec


COMP_CASES = [
    ('water_percent', {'H': 11.19, 'O': 88.81}, 'percent', None),
    ('glucose_percent', {'C': 40.00, 'H': 6.71, 'O': 53.29}, 'percent', 180.156),
    ('benzene_percent', {'C': 92.26, 'H': 7.74}, 'percent', 78.11),
    ('ethanol_percent', {'C': 52.14, 'H': 13.13, 'O': 34.73}, 'percent', 46.068),
    ('iron_oxide_masses', {'Fe': 111.69, 'O': 47.997}, 'masses', None),
    ('acetic_acid_masses', {'C': 24.022, 'H': 4.032, 'O': 31.998}, 'masses', 60.052),
    ('noisy_ratio', {'C': 1.0, 'H': 1.0}, 'masses', None),
    ('h_only', {'H': 1.008}, 'masses', None),
]
for slug, comp, mode, mm in COMP_CASES:
    r = empirical(comp, mode, mm)
    print(f'comp.{slug}.formula\t{r["formula"]}\tk={r["kmult"]}\tdev={r["dev"]!r}\texact={r["exact"]}')
    print(f'comp.{slug}.mass\t{r["mass"]!r}')
    for el in sorted(r['atoms']):
        print(f'comp.{slug}.atoms.{el}\tcount={r["atoms"][el]}\tmoles={r["moles"][el]!r}'
              f'\tratio={r["ratios"][el]!r}\tpct={r["pct"][el]!r}')
    if 'molformula' in r:
        print(f'comp.{slug}.molecular\t{r["molformula"]}\tn={r["nearest"]}\tnrel={r["nrel"]!r}'
              f'\tnmoles={r["nmol"]!r}\texact={r["nexact"]}\timplied={r["implied"]!r}\tmerr={r["merr"]!r}')

# ---------------------------------------------------------------------------
# equilibrium anchors (weak acid / weak base ICE, exact quadratic)
# Independent of the TS implementation: the positive root of x^2 + K x - K C = 0
# written out by hand, pH from the definition. Published textbook values are the
# ruler (0.1 M acetic acid pH 2.87, 0.1 M ammonia pH 11.13).
print()
print('# --- equilibrium anchors (weak acid / weak base ICE) ---')
KW = 1.0e-14


def weak_x(k, c):
    return (-k + math.sqrt(k * k + 4.0 * k * c)) / 2.0


def weak_case(kind, k, c):
    x = weak_x(k, c)
    pK = -math.log10(k)
    if kind == 'weak_acid':
        h, oh = x, KW / x
        pH = -math.log10(h)
        pOH = -math.log10(oh)
    else:
        oh, h = x, KW / x
        pOH = -math.log10(oh)
        pH = -math.log10(h)
    alpha = x / c
    approx = math.sqrt(k * c)
    return dict(kind=kind, k=k, c=c, x=x, pK=pK, h=h, oh=oh, pH=pH, pOH=pOH, alpha=alpha,
                approx=approx, approx_err=(approx - x) / x * 100.0,
                k_back=x * x / (c - x), oh_from_kw=KW / h)


EQUIL_CASES = [
    ('acetic_0.1', 'weak_acid', 1.8e-5, 0.1),
    ('acetic_0.01', 'weak_acid', 1.8e-5, 0.01),
    ('acetic_1.0', 'weak_acid', 1.8e-5, 1.0),
    ('formic_0.05', 'weak_acid', 1.8e-4, 0.05),
    ('ammonia_0.1', 'weak_base', 1.8e-5, 0.1),
    ('ammonia_0.001', 'weak_base', 1.8e-5, 0.001),
]
for slug, kind, k, c in EQUIL_CASES:
    r = weak_case(kind, k, c)
    print(f'eq.{slug}\tkind={kind}\tk={r["k"]!r}\tc={r["c"]!r}\tx={r["x"]!r}\tpK={r["pK"]!r}')
    print(f'eq.{slug}.conc\tH={r["h"]!r}\tOH={r["oh"]!r}\tpH={r["pH"]!r}\tpOH={r["pOH"]!r}')
    print(f'eq.{slug}.alpha\talpha={r["alpha"]!r}\tpercent={r["alpha"] * 100!r}')
    print(f'eq.{slug}.approx\tx={r["approx"]!r}\terr={r["approx_err"]!r}')
    print(f'eq.{slug}.check\tk_back={r["k_back"]!r}\toh_from_kw={r["oh_from_kw"]!r}')

# ---------------------------------------------------------------------------
# ksp anchors (Ksp -> molar solubility, ion grouping, common-ion effect)
# Independent of the TS implementation: the salt is matched against the ion
# table by ELEMENT COUNTS (a*cat + b*an == salt) instead of splitting the
# parsed group list, and the common-ion case is solved by Newton iteration on
# the residual instead of bisection. Published handbook values are the ruler
# (AgCl 1.8e-10 -> 1.3e-5 M, Ca(OH)2 5.5e-6 -> 1.1e-2 M, PbCl2 1.7e-5 ->
# 1.6e-2 M, Ag2CrO4 1.1e-12 -> 6.5e-5 M).
print()
print('# --- ksp anchors (solubility product -> molar solubility) ---')

KS_CATIONS = [
    ('Li^+', 1, 'Li'), ('Na^+', 1, 'Na'), ('K^+', 1, 'K'), ('Rb^+', 1, 'Rb'), ('Cs^+', 1, 'Cs'),
    ('Ag^+', 1, 'Ag'), ('Cu^+', 1, 'Cu'), ('Hg^+', 1, 'Hg'), ('NH4^+', 1, 'NH4'),
    ('Be^2+', 2, 'Be'), ('Mg^2+', 2, 'Mg'), ('Ca^2+', 2, 'Ca'), ('Sr^2+', 2, 'Sr'), ('Ba^2+', 2, 'Ba'),
    ('Zn^2+', 2, 'Zn'), ('Cd^2+', 2, 'Cd'), ('Hg^2+', 2, 'Hg'), ('Hg2^2+', 2, 'Hg2'),
    ('Cu^2+', 2, 'Cu'), ('Fe^2+', 2, 'Fe'), ('Co^2+', 2, 'Co'), ('Ni^2+', 2, 'Ni'), ('Mn^2+', 2, 'Mn'),
    ('Pb^2+', 2, 'Pb'), ('Sn^2+', 2, 'Sn'),
    ('Al^3+', 3, 'Al'), ('Ga^3+', 3, 'Ga'), ('In^3+', 3, 'In'), ('Cr^3+', 3, 'Cr'), ('Fe^3+', 3, 'Fe'),
    ('Bi^3+', 3, 'Bi'),
    ('Sn^4+', 4, 'Sn'),
]
KS_ANIONS = [
    ('F^-', -1, 'F'), ('Cl^-', -1, 'Cl'), ('Br^-', -1, 'Br'), ('I^-', -1, 'I'), ('OH^-', -1, 'OH'),
    ('NO3^-', -1, 'NO3'), ('NO2^-', -1, 'NO2'), ('HCO3^-', -1, 'CHO3'), ('HSO4^-', -1, 'HSO4'),
    ('H2PO4^-', -1, 'H2PO4'), ('ClO^-', -1, 'ClO'), ('ClO3^-', -1, 'ClO3'), ('ClO4^-', -1, 'ClO4'),
    ('MnO4^-', -1, 'MnO4'), ('CN^-', -1, 'CN'), ('SCN^-', -1, 'CNS'), ('CH3COO^-', -1, 'C2H3O2'),
    ('IO3^-', -1, 'IO3'), ('BrO3^-', -1, 'BrO3'),
    ('O^2-', -2, 'O'), ('S^2-', -2, 'S'), ('Se^2-', -2, 'Se'), ('O2^2-', -2, 'O2'),
    ('SO4^2-', -2, 'SO4'), ('SO3^2-', -2, 'SO3'), ('S2O3^2-', -2, 'S2O3'), ('CO3^2-', -2, 'CO3'),
    ('HPO4^2-', -2, 'HPO4'), ('C2O4^2-', -2, 'C2O4'), ('CrO4^2-', -2, 'CrO4'), ('Cr2O7^2-', -2, 'Cr2O7'),
    ('SiO3^2-', -2, 'SiO3'),
    ('PO4^3-', -3, 'PO4'), ('PO3^3-', -3, 'PO3'), ('AsO4^3-', -3, 'AsO4'),
]


def _ion_counts(formula):
    _m, counts, _h, charge, _c = parse_formula(formula)
    if charge is not None:
        raise ValueError(f'ion formula {formula} must be neutral')
    return counts


KS_IONS = ([('cat', label, q, _ion_counts(f)) for label, q, f in KS_CATIONS]
           + [('an', label, q, _ion_counts(f)) for label, q, f in KS_ANIONS])


def _dicts_add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
    return out


def _scaled(counts, k):
    return {s: v * k for s, v in counts.items()}


def ksp_groupings(salt):
    """Every (cation, count, anion, count) whose ions reconstruct the salt exactly."""
    _m, salt_counts, hydrate, charge, coeff = parse_formula(salt)
    if hydrate is not None or charge is not None or coeff is not None:
        raise ValueError(f'{salt} is not a plain neutral formula')
    found = []
    for kind_c, label_c, qc, comp_c in KS_IONS:
        if kind_c != 'cat':
            continue
        for kind_a, label_a, qa, comp_a in KS_IONS:
            if kind_a != 'an':
                continue
            for a in range(1, 13):
                # electroneutrality fixes b for each a
                num = a * qc
                if num % (-qa) != 0:
                    continue
                b = num // (-qa)
                if b < 1 or b > 12:
                    continue
                if _dicts_add(_scaled(comp_c, a), _scaled(comp_a, b)) == salt_counts:
                    found.append((label_c, a, label_a, b, comp_c, comp_a))
    return found, salt_counts


def ksp_case(salt, ksp, common=None, common_c=None, force_cation=None, force_anion=None) -> dict:
    found, salt_counts = ksp_groupings(salt)
    if force_cation is not None or force_anion is not None:
        found = [f for f in found
                 if (force_cation is None or f[0] == force_cation)
                 and (force_anion is None or f[2] == force_anion)]
    if len(found) != 1:
        raise ValueError(f'{salt}: expected exactly one grouping, got {[f[:4] for f in found]}')
    label_c, a, label_a, b, _cc, _ca = found[0]
    qc = [q for k, l, q, _ in KS_IONS if k == 'cat' and l == label_c][0]
    qa = [q for k, l, q, _ in KS_IONS if k == 'an' and l == label_a][0]
    factor = (a ** a) * (b ** b)
    order = a + b
    s = (ksp / factor) ** (1.0 / order)
    mass = parse_formula(salt)[0]
    out = dict(salt=salt, ksp=ksp, cation=label_c, a=a, anion=label_a, b=b,
               charge_balance=a * qc + b * qa, factor=factor, order=order,
               s=s, mass_solubility=s * mass, molar_mass=mass,
               cat_conc=a * s, an_conc=b * s, common=None)
    if common is not None:
        # Newton on f(s) = (cat-side)(an-side) - Ksp, starting from Ksp/C.
        cation_is_common = (common == label_c)
        present = common_c

        def f(x):
            cat = a * x + (present if cation_is_common else 0.0)
            an = b * x + (0.0 if cation_is_common else present)
            return (cat ** a) * (an ** b) - ksp

        x = max(ksp / present, 1e-300)
        for _ in range(200):
            h = x * 1e-7 + 1e-30
            d = (f(x + h) - f(x - h)) / (2 * h)
            if d == 0:
                break
            step = f(x) / d
            x -= step
            if x < 0:
                x = 0.0
            if abs(step) <= 1e-16 * max(x, 1e-300):
                break
        out['common'] = dict(symbol=common, c=present, s=x, reduction=s / x,
                            residual=f(x) / ksp, mass_solubility=x * mass)
    return out


KSP_CASES = [
    ('AgCl', 1.8e-10, None, None),
    ('Ca(OH)2', 5.5e-6, None, None),
    ('PbCl2', 1.7e-5, None, None),
    ('Ag2CrO4', 1.1e-12, None, None),
    ('BaSO4', 1.1e-10, None, None),
    ('Mg(OH)2', 5.6e-12, None, None),
    ('Ca3(PO4)2', 2.07e-33, None, None),
    ('Fe(OH)3', 2.79e-39, None, None),
    ('LiF', 1.8e-3, None, None),
    ('NH4NO3', 1.0e-3, None, None),
    # mercury(I) is written Hg2^2+ in solution chemistry and Hg+ in elementary
    # texts; both readings reconstruct Hg2Cl2, so the tool must be told which.
    ('Hg2Cl2', 1.3e-18, None, None, 'Hg2^2+', None),
    ('AgCl', 1.8e-10, 'Cl^-', 0.01),
    ('AgCl', 1.8e-10, 'Ag^+', 0.01),
]
for case in KSP_CASES:
    salt, ksp, common, common_c = case[:4]
    force_cation, force_anion = (case[4], case[5]) if len(case) > 4 else (None, None)
    try:
        r = ksp_case(salt, ksp, common, common_c, force_cation, force_anion)
    except ValueError as err:
        print(f'ksp.{salt}.{ksp!r}.{common}\tERROR {err}')
        continue
    slug = f'{salt}_{ksp!r}_{common}'.replace('^', '').replace('(', '').replace(')', '').replace(' ', '')
    print(f'ksp.{slug}\tcations={r["cation"]}x{r["a"]}\tanions={r["anion"]}x{r["b"]}'
          f'\tcharge_balance={r["charge_balance"]}\tfactor={r["factor"]!r}\torder={r["order"]}')
    print(f'ksp.{slug}.s\tmolar_solubility={r["s"]!r}\tmass_solubility={r["mass_solubility"]!r}'
          f'\tmolar_mass={r["molar_mass"]!r}')
    print(f'ksp.{slug}.conc\tcat={r["cat_conc"]!r}\tan={r["an_conc"]!r}')
    if r['common']:
        c = r['common']
        print(f'ksp.{slug}.common\ts={c["s"]!r}\treduction={c["reduction"]!r}\tresidual={c["residual"]!r}'
              f'\tmass_solubility={c["mass_solubility"]!r}')

# ---------------------------------------------------------------------------
# --check: self-verification. The ruler must be OUTSIDE this file's own maths
# (published textbook values + definitional identities) — an oracle that only
# agrees with the implementation would pass even when both are wrong (R63).
def self_check():
    failures = []

    def check(label, got, want, tol):
        if got is None or not (abs(got - want) <= tol):
            failures.append(f'{label}: {got!r} != {want!r} (tol {tol})')

    for formula, published in [('H2O', 18.015), ('H2SO4', 98.08), ('NaCl', 58.44), ('Ca(OH)2', 74.09),
                               ('CuSO4·5H2O', 249.68), ('C6H12O6', 180.16), ('KMnO4', 158.03), ('NH3', 17.03)]:
        check(f'molar mass {formula}', _fmass(formula), published, 0.02)

    for eq, want in [('C3H8 + O2 -> CO2 + H2O', [1, 5, 3, 4]),
                     ('KMnO4 + HCl -> KCl + MnCl2 + Cl2 + H2O', [2, 16, 2, 2, 5, 8])]:
        kind, coeffs, _u, _e = balance_equation(eq)
        if kind != 'ok' or coeffs != want:
            failures.append(f'balance {eq}: {kind} {coeffs} != {want}')

    check('R in L·atm/(mol·K)', GAS_R, 0.082057366, 1e-8)
    check('molar volume at 0 °C, 1 atm', GAS_R * 273.15, 22.41396954, 1e-6)
    check('pH of 0.01 M strong acid', -math.log10(0.01), 2.0, 1e-12)
    check('pH of 0.01 M strong base', 14 + math.log10(0.01), 12.0, 1e-12)

    for slug, eq, amounts, product, actual in STOICH_CASES:
        extent, reactants, products, limiting, _y = stoich(eq, amounts, product, actual)
        if extent <= 0:
            failures.append(f'stoich {slug}: non-positive extent {extent!r}')
        mass_in = sum(r[5] * r[2] for r in reactants)
        mass_out = sum(p[3] * p[2] for p in products)
        if abs(mass_in - mass_out) > 1e-9 * max(1.0, mass_in):
            failures.append(f'stoich {slug}: mass not conserved ({mass_in!r} vs {mass_out!r})')
        for r in reactants:
            if r[6] < -1e-12:
                failures.append(f'stoich {slug}: negative leftover for {r[0]}')
            if r[5] > r[3] + 1e-12:
                failures.append(f'stoich {slug}: consumed more {r[0]} than available')
        if not limiting:
            failures.append(f'stoich {slug}: no limiting reagent')

    # published: exactly 1 mol propane : 5 mol O2 is the stoichiometric ratio —
    # both reactants run out together and nothing is left over.
    extent, reactants, products, limiting, _y = stoich(
        'C3H8 + O2 -> CO2 + H2O', {'C3H8': {'moles': 1.0}, 'O2': {'moles': 5.0}})
    if sorted(limiting) != ['C3H8', 'O2']:
        failures.append(f'propane exact ratio: limiting {limiting}')
    if any(r[6] > 1e-12 for r in reactants):
        failures.append('propane exact ratio: leftover present')
    # published: CO2 = 44.01 g/mol → 3 mol CO2 = 132.03 g
    co2 = [p for p in products if p[0] == 'CO2'][0]
    check('3 mol CO2 mass', co2[4], 132.03, 0.02)
    # published: 1 mol Zn (65.38 g/mol) is 65.38 g and liberates 1 mol H2 (2.016 g)
    _e, r_zn, p_h2, _l, _y = stoich('Zn + 2HCl -> ZnCl2 + H2',
                                    {'Zn': {'grams': 65.38}, 'HCl': {'moles': 2.0}})
    check('1 mol Zn mass', r_zn[0][3] * r_zn[0][2], 65.38, 0.02)
    check('1 mol H2 mass', [p for p in p_h2 if p[0] == 'H2'][0][4], 2.016, 0.01)

    # --- composition: published formulas are the ruler (percent data from the
    # textbook analyses of water / glucose / benzene / ethanol, masses for the
    # oxides and acetic acid) ---
    for slug, comp, mode, mm, want_formula, want_mol in [
        ('water', {'H': 11.19, 'O': 88.81}, 'percent', None, 'H2O', None),
        ('glucose', {'C': 40.00, 'H': 6.71, 'O': 53.29}, 'percent', 180.156, 'CH2O', 'C6H12O6'),
        ('benzene', {'C': 92.26, 'H': 7.74}, 'percent', 78.11, 'CH', 'C6H6'),
        ('ethanol', {'C': 52.14, 'H': 13.13, 'O': 34.73}, 'percent', 46.068, 'C2H6O', 'C2H6O'),
        ('iron_oxide', {'Fe': 111.69, 'O': 47.997}, 'masses', None, 'Fe2O3', None),
        ('acetic_acid', {'C': 24.022, 'H': 4.032, 'O': 31.998}, 'masses', 60.052, 'CH2O', 'C2H4O2'),
    ]:
        r = empirical(comp, mode, mm)
        if r['formula'] != want_formula or not r['exact']:
            failures.append(f'composition {slug}: {r["formula"]} (exact={r["exact"]}) != {want_formula}')
        if want_mol is not None and r.get('molformula') != want_mol:
            failures.append(f'composition {slug}: molecular {r.get("molformula")} != {want_mol}')
        # invariant: the round trip must reproduce the input percentages
        if mode == 'percent':
            for el, given in comp.items():
                check(f'composition {slug} pct {el}', r['pct'][el], given, 0.5)
        # invariant: empirical subscripts are coprime by construction
        g = 0
        for n in r['atoms'].values():
            g = math.gcd(g, n)
        if g != 1:
            failures.append(f'composition {slug}: subscripts not coprime (gcd {g})')
        if want_mol is not None and abs(r['nearest'] * r['mass'] - mm) / mm > 0.02:
            failures.append(f'composition {slug}: molecular multiplier off by {r["merr"]:.3f} %')
    noisy = empirical({'C': 1.0, 'H': 1.0}, 'masses')
    if noisy['exact']:
        failures.append('composition noisy_ratio: must be flagged inexact (dev 0.084 > 0.05)')

    # --- equilibrium: definitional identity x^2/(C-x) = K and published pH ---
    for kind, k, c, want_pH in [
        ('weak_acid', 1.8e-5, 0.1, 2.87),
        ('weak_acid', 1.8e-4, 0.05, 2.54),
        ('weak_base', 1.8e-5, 0.1, 11.13),
    ]:
        r = weak_case(kind, k, c)
        if abs(r['k_back'] - k) > 1e-9 * k:
            failures.append(f'equilibrium {kind} {c}: x^2/(C-x) = {r["k_back"]!r} != K {k!r}')
        check(f'equilibrium pH {kind} {c} M', r['pH'], want_pH, 0.01)
        if not (0 < r['x'] < c):
            failures.append(f'equilibrium {kind} {c}: x={r["x"]!r} outside (0, C)')
    errs = [(c, weak_case('weak_acid', 1.8e-5, c)['approx_err']) for c in (0.01, 0.1, 1.0)]
    for (c1, e1), (c2, e2) in zip(errs, errs[1:]):
        if not e2 < e1:
            failures.append(f'equilibrium: sqrt(KC) error not decreasing ({c1}: {e1:.3f} %, {c2}: {e2:.3f} %)')

    # --- ksp: published molar solubilities are the ruler, plus the structural
    # identities (charge balance, Ksp reconstruction, common ion shifts it down)
    for salt, ksp, want_s in [
        ('AgCl', 1.8e-10, 1.3e-5),
        ('Ca(OH)2', 5.5e-6, 1.1e-2),
        ('PbCl2', 1.7e-5, 1.6e-2),
        ('Ag2CrO4', 1.1e-12, 6.5e-5),
        ('BaSO4', 1.1e-10, 1.0e-5),
        ('Mg(OH)2', 5.6e-12, 1.1e-4),
    ]:
        r = ksp_case(salt, ksp)
        check(f'ksp {salt} molar solubility', r['s'], want_s, want_s * 0.05)
        if r['charge_balance'] != 0:
            failures.append(f'ksp {salt}: net charge {r["charge_balance"]} != 0')
        if abs(r['factor'] * r['s'] ** r['order'] - ksp) > 1e-12 * ksp:
            failures.append(f'ksp {salt}: factor*s^N = {r["factor"] * r["s"] ** r["order"]!r} != Ksp {ksp!r}')
    # Hg2Cl2 must be ambiguous BEFORE an ion is forced: both readings reconstruct it
    hg_all, _ = ksp_groupings('Hg2Cl2')
    hg_labels = sorted({(f[0], f[2]) for f in hg_all})
    if hg_labels != [('Hg2^2+', 'Cl^-'), ('Hg^+', 'Cl^-')]:
        failures.append(f'ksp Hg2Cl2: expected two readings Hg2^2+ / Hg^+, got {hg_labels}')
    # common ion: exact identity s*(a s + C)^... = Ksp and a strictly smaller s
    for common, c in [('Cl^-', 0.01), ('Ag^+', 0.01)]:
        r = ksp_case('AgCl', 1.8e-10, common, c)
        cm = r['common']
        if abs(cm['s'] * (cm['s'] + c) - 1.8e-10) > 1e-12 * 1.8e-10:
            failures.append(f'ksp AgCl + {common}: (s)(s+C) = {cm["s"] * (cm["s"] + c)!r} != Ksp')
        if not cm['s'] < r['s']:
            failures.append(f'ksp AgCl + {common}: common ion must lower solubility')
        if abs(cm['s'] - 1.8e-8) > 0.5e-8:
            failures.append(f'ksp AgCl + {common}: s = {cm["s"]!r}, published estimate 1.8e-8')
    # monotonic in Ksp, and the exponent really is a+b
    lo = ksp_case('Ca3(PO4)2', 2.07e-33)['s']
    hi = ksp_case('Ca3(PO4)2', 2.07e-30)['s']
    if not hi > lo:
        failures.append(f'ksp monotonicity: s did not rise with Ksp ({lo!r} -> {hi!r})')
    if abs((hi / lo) - 10 ** (3.0 / 5.0)) > 1e-9 * 10 ** (3.0 / 5.0):
        failures.append(f'ksp Ca3(PO4)2: 1000x Ksp must raise s by 10^(3/5), got {hi / lo!r}')

    if failures:
        print('SELF-CHECK FAILED')
        for line in failures:
            print('  [x]', line)
        return 1
    print('SELF-CHECK OK: 8 published molar masses, 2 textbook balances, R/molar volume, '
          'pH definition, 7 stoichiometry cases (mass conservation + leftover bounds), '
          '6 published compositions (empirical + molecular formulas), '
          '3 published weak acid/base pH values + x^2/(C-x) == K, '
          '6 published Ksp solubility values + charge balance/Ksp reconstruction + '
          'Hg2Cl2 two readings + common-ion identity')
    return 0


if '--check' in sys.argv:
    sys.exit(self_check())
