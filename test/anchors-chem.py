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
import math, re

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
