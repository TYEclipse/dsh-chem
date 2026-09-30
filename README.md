# dsh-chem

Deterministic chemistry math for [DeepSeek Harness](https://github.com/deepseek-ai/deepseek-harness) (`dsh`) — **zero runtime dependencies**, pure arithmetic over IUPAC data.

> 中文简介：化学计量数学工具箱。六个零依赖确定性工具——化学式解析与摩尔质量（支持嵌套括号、结晶水、离子电荷）、质量/摩尔/粒子数换算、溶液浓度（摩尔浓度、稀释 C1V1=C2V2、百分比浓度）、pH/pOH 与缓冲溶液、配平化学方程式（精确有理数运算，含电荷守恒）、理想气体定律（PV=nRT 求解、摩尔体积、气体密度与由密度反推摩尔质量）。专门替代 LLM 高频出错的化学心算。

## Why this plugin

Agents get chemistry arithmetic wrong in predictable ways:

- **Molar-mass sums** — miscounted subscripts, dropped hydrates (CuSO4·5H2O is 249.68 g/mol, not 159.61).
- **Element case sensitivity** — "Co" is cobalt (58.93), "CO" is carbon monoxide (28.01).
- **Dilution math** — C1V1 = C2V2 with a wrong or missing factor.
- **pH log/antilog arithmetic** — pH of 0.01 M HCl is 2, not 12.

These tools replace all of it with deterministic tables and formulas. No network, no state, no dependencies — install and call.

## Tools

### `formula_mass` — parse a formula and compute molar mass

Accepts nested parentheses, hydrate dots (`CuSO4·5H2O`, `*` also works), leading coefficients (`2H2O`) and ionic charges. Uses IUPAC standard atomic weights; elements without stable isotopes use the most-stable-isotope mass number (flagged with `usedMassNumber`).

```json
{ "formula": "H2SO4" }
```

```json
{
  "valid": true,
  "formula": "H2SO4",
  "molarMass": 98.072,
  "elements": [
    { "symbol": "H", "name": "Hydrogen", "count": 2, "atomicMass": 1.008, "massContribution": 2.016, "massFraction": 0.02056 },
    { "symbol": "O", "name": "Oxygen", "count": 4, "atomicMass": 15.999, "massContribution": 63.996, "massFraction": 0.65254 },
    { "symbol": "S", "name": "Sulfur", "count": 1, "atomicMass": 32.06, "massContribution": 32.06, "massFraction": 0.32690 }
  ],
  "usedMassNumber": false
}
```

**Charge disambiguation rule** (documented, deterministic): digits immediately before a trailing `+`/`-` belong to the *ionic charge* only when every element parsed so far is the same single element — `Fe3+` is Fe with 3+, `S2-` is S with 2-, `Ca2+` is Ca with 2+. Otherwise they are a subscript and the bare sign is a 1-unit charge — `NH4+` is N + H4 with `+`, `NO3-` is N + O3 with `-`. The explicit form `SO4^2-` is always the charge, and `O2^+` escapes the rule for subscript-before-charge. Single lowercase symbols are auto-corrected (`h2o` → H2O); multi-letter lowercase runs are rejected with a case-sensitivity hint.

### `mass_moles` — mass / amount / particle conversions

Six modes (`mass_to_moles`, `moles_to_mass`, `mass_to_particles`, `particles_to_mass`, `moles_to_particles`, `particles_to_moles`), molar mass from a formula or given directly, Avogadro constant 6.02214076e23.

```json
{ "mode": "mass_to_moles", "value": 98, "formula": "H2SO4" }
```

```json
{ "valid": true, "mode": "mass_to_moles", "molarMass": 98.072, "avogadro": 6.02214076e+23, "formula": "H2SO4", "moles": 0.9992658455012644 }
```

### `solution` — concentration math

Four operations:

- `molarity` — mass of solute + volume → mol/L: `{ "operation": "molarity", "massG": 58.44, "volumeMl": 250, "formula": "NaCl" }` → `molarity: 4`
- `prepare` — target molarity + volume → solute mass: `{ "operation": "prepare", "molarity": 0.1, "volumeMl": 500, "formula": "NaOH" }` → `massG: 1.99985`
- `dilute` — C1V1 = C2V2, give any three of `c1`/`v1`/`c2`/`v2` (M, mL) and the missing one is computed: `{ "operation": "dilute", "c1": 2, "c2": 0.5, "v2": 250 }` → `v1: 62.5`, `solved: "v1"`
- `percent` — `wv` (g per 100 mL), `vv` (mL per 100 mL) or `mass` percent: `{ "operation": "percent", "percentType": "wv", "soluteG": 5, "solutionMl": 200 }` → `percent: 2.5`

### `ph` — pH chemistry at 25 °C (Kw = 1e-14)

Five operations:

- `from_h` — `{ "operation": "from_h", "hConc": 0.001 }` → `pH: 3, pOH: 11, ohConc: 1e-11`
- `from_oh` — [OH-] → pH/pOH
- `from_ph` — pH → [H+]
- `from_poh` — pOH → [OH-]
- `buffer` — Henderson–Hasselbalch: `{ "operation": "buffer", "pKa": 4.76, "acidConc": 0.2, "baseConc": 0.1 }` → `pH: 4.4589700043360185, bufferRatio: 0.5`

Negative pH for >1 M acids is supported ([H+] = 2 → pH −0.301).

### `balance_equation` — balance an equation (exact integer coefficients)

The solver runs entirely over BigInt-backed rationals, so coefficients are exact integers — never `2.9999997`. Write species separated by ` + ` and one arrow (`->`, `=`, `→`); input coefficients like `2H2` are accepted and ignored (scaling a column cannot change the solution).

```json
{ "equation": "C3H8 + O2 -> CO2 + H2O" }
```

```json
{
  "valid": true,
  "balanced": true,
  "unique": true,
  "equation": "C3H8 + 5 O2 -> 3 CO2 + 4 H2O",
  "reactants": [{ "formula": "C3H8", "coefficient": 1 }, { "formula": "O2", "coefficient": 5 }],
  "products": [{ "formula": "CO2", "coefficient": 3 }, { "formula": "H2O", "coefficient": 4 }],
  "elements": [
    { "symbol": "C", "left": 3, "right": 3 },
    { "symbol": "H", "left": 8, "right": 8 },
    { "symbol": "O", "left": 10, "right": 10 }
  ],
  "chargeBalanced": false
}
```

Deterministic behaviour worth knowing:

- **Charged species force charge conservation too** — `MnO4- + Fe2+ + H+ -> Mn2+ + Fe3+ + H2O` → 1/5/8/1/5/4 with `chargeBalanced: true`.
- **Hydrate water counts as water** — `CuSO4·5H2O -> CuSO4 + H2O` → 1/1/5.
- **Separation rule**: a bare `+` next to a symbol is a *charge*, so `Na+ + Cl-` keeps its ions while the compact form `Fe+O2->Fe2O3` splits. Charged species must use the spaced form (the compact form is ambiguous by construction).
- **Honest failures**: `H2 + O2 -> H2O + N2` returns `valid: true, balanced: false` with the reason (“at least one species would need coefficient 0”) instead of silently dropping N2. Equations with several independent balances (`Fe + O2 -> Fe2O3 + Fe3O4`) return `unique: false` plus the smallest solution by total coefficient sum. Malformed input (no arrow, unknown element) returns `valid: false`.

### `gas_law` — ideal gas arithmetic with unit conversions

Fixed units, no unit guessing: pressure **atm** (also reported in kPa/mmHg), volume **L** (also mL), amount **mol**, temperature **K** (also °C). Four operations:

- `ideal` — give exactly three of `pressureAtm`, `volumeL`, `moles`, `temperatureK` and the fourth is solved from PV = nRT: `{ "operation": "ideal", "moles": 1, "temperatureK": 273.15, "volumeL": 22.4 }` → `pressureAtm: 1.0006236404024167`.
- `molar_volume` — V/n = RT/P, defaults 0 °C / 1 atm: `{ "operation": "molar_volume" }` → `molarVolumeLPerMol: 22.413969545014133` (textbook 22.414).
- `density` — ρ = P·M/(R·T) from `formula` or `molarMass`: `{ "operation": "density", "formula": "CO2", "pressureAtm": 1, "temperatureK": 273.15 }` → `densityGPerL: 1.9634630051413435`.
- `molar_mass` — M = ρ·R·T/P from a measured density: `{ "operation": "molar_mass", "densityGPerL": 1.96, "pressureAtm": 1, "temperatureK": 273.15 }` → `molarMass: 43.9313803082277`.

R = 8.31446261815324 / 101.325 L·atm/(mol·K). Missing operands and non-positive values are reported as errors rather than guessed.

## Install

```bash
dsh plugin --profile <profile> add github:TYEclipse/dsh-chem
```

(Requires `pnpm` on PATH — `dsh plugin` forwards installs to it.)

## Verify

```bash
dsh plugin --profile <profile> list      # dsh-chem present
dsh --dump-config --profile <profile>    # == dsh-chem layer present
```

## Development

```bash
pnpm install
pnpm build     # tsc -> dist/ (committed; git installs do not build)
pnpm test      # vitest, 136 tests
pnpm lint      # oxlint src test
```

All numeric anchors are generated by the independent Python oracle `test/anchors-chem.py` (own periodic table + own tokenizer + own Fraction-based equation balancer) and cross-validated against published textbook values (H2O 18.015, H2SO4 98.08, NaCl 58.44, Ca(OH)2 74.09, CuSO4·5H2O 249.68, C6H12O6 180.16, KMnO4 158.03, NH3 17.03; propane 1/5/3/4; permanganate/HCl 2/16/2/2/5/8; 22.414 L/mol at STP). Run it with `python3 test/anchors-chem.py`.

## License

MIT
