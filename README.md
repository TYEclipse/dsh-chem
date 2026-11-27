# dsh-chem

Deterministic chemistry math for [DeepSeek Harness](https://github.com/deepseek-ai/deepseek-harness) (`dsh`) — **zero runtime dependencies**, pure arithmetic over IUPAC data.

> 中文简介：化学计量数学工具箱。十个零依赖确定性工具——化学式解析与摩尔质量（支持嵌套括号、结晶水、离子电荷）、质量/摩尔/粒子数换算、溶液浓度（摩尔浓度、稀释 C1V1=C2V2、百分比浓度）、pH/pOH 与缓冲溶液、配平化学方程式（精确有理数运算，含电荷守恒）、理想气体定律（PV=nRT 求解、摩尔体积、气体密度与由密度反推摩尔质量）、限量试剂与产率、元素分析反推化学式、弱酸弱碱 ICE 平衡、溶度积 Ksp 求摩尔溶解度（含同离子效应）。专门替代 LLM 高频出错的化学心算。

## Why this plugin

Agents get chemistry arithmetic wrong in predictable ways:

- **Molar-mass sums** — miscounted subscripts, dropped hydrates (CuSO4·5H2O is 249.68 g/mol, not 159.61).
- **Element case sensitivity** — "Co" is cobalt (58.93), "CO" is carbon monoxide (28.01).
- **Dilution math** — C1V1 = C2V2 with a wrong or missing factor.
- **pH log/antilog arithmetic** — pH of 0.01 M HCl is 2, not 12.
- **Limiting-reagent reasoning** — dividing two amounts by the coefficients and comparing, then carrying the extent of reaction through every product.

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

### `stoichiometry` — limiting reagent, yield and leftovers

Give the equation plus how much of **each** reactant you actually have (moles or grams) and the tool
answers the questions that follow: the extent of reaction, which reactant runs out first, how much of
every product forms, what is left over — and, with `product` + `actualYieldG`, the theoretical and
percent yield.

```json
{ "equation": "C3H8 + O2 -> CO2 + H2O",
  "reactants": [{ "formula": "C3H8", "moles": 0.5 }, { "formula": "O2", "moles": 2 }],
  "product": "CO2", "actualYieldG": 44 }
```

```json
{
  "valid": true,
  "equation": "C3H8 + 5O2 -> 3CO2 + 4H2O",
  "extent": 0.4,
  "limiting": ["O2"],
  "reactants": [
    { "formula": "C3H8", "coefficient": 1, "molarMass": 44.097, "availableMoles": 0.5,
      "availableGrams": 22.0485, "consumedMoles": 0.4, "leftoverMoles": 0.1,
      "leftoverGrams": 4.4097, "limiting": false },
    { "formula": "O2", "coefficient": 5, "molarMass": 31.998, "availableMoles": 2,
      "availableGrams": 63.996, "consumedMoles": 2, "leftoverMoles": 0,
      "leftoverGrams": 0, "limiting": true }
  ],
  "products": [
    { "formula": "CO2", "coefficient": 3, "molarMass": 44.009, "moles": 1.2, "grams": 52.8108 },
    { "formula": "H2O", "coefficient": 4, "molarMass": 18.015, "moles": 1.6, "grams": 28.824 }
  ],
  "yieldProduct": "CO2",
  "theoreticalYieldMol": 1.2,
  "theoreticalYieldG": 52.8108,
  "actualYieldG": 44,
  "percentYield": 83.3162913646
}
```

Deterministic behaviour worth knowing:

- **Every reactant needs an amount** — an unstated amount would make the limiting reagent undefined, so it is an error naming the missing species instead of an assumption about “excess”.
- **Co-limiting is reported, not broken by a tie-break**: at exactly 2 mol H2 : 1 mol O2 both are returned in `limiting`, and nothing is left over.
- **Grams are converted with the same IUPAC molar masses** `formula_mass` reports, so `6.54 g Zn` is the same number everywhere.
- **Coefficients come from the exact rational balancer** — this tool balances first (`C3H8 + O2 -> …` works as written), charge and hydrate included.
- **Honest failures**: a species that is not a reactant, both/neither of `moles`/`grams`, a duplicate entry, a non-positive amount, a product that is not in the equation, or `actualYieldG` without `product` all return `valid: false` with the reason.

### `formula_from_composition` — elemental analysis → empirical (and molecular) formula

Give the mass percent (or the grams) of each element and get the empirical formula, the integer
multiplier the analysis needed, and — when `molarMass` is also given — the molecular formula.

```json
{ "mode": "percent",
  "parts": [{ "symbol": "C", "value": 40.0 }, { "symbol": "H", "value": 6.71 }, { "symbol": "O", "value": 53.29 }],
  "molarMass": 180.156 }
```

```json
{
  "valid": true,
  "mode": "percent",
  "empiricalFormula": "CH2O",
  "empiricalMolarMass": 30.026,
  "multiplier": 1,
  "maxDeviation": 0.00114558531746,
  "exact": true,
  "elements": [
    { "symbol": "C", "massG": 40, "atomicMass": 12.011, "moles": 3.33028057614, "ratio": 1, "subscript": 1, "massFraction": 40.0019982682 },
    { "symbol": "H", "massG": 6.71, "atomicMass": 1.008, "moles": 6.65674603175, "ratio": 1.99885441468, "subscript": 2, "massFraction": 6.7141810431 },
    { "symbol": "O", "massG": 53.29, "atomicMass": 15.999, "moles": 3.33083317707, "ratio": 1.00016593225, "subscript": 1, "massFraction": 53.2838206887 }
  ],
  "molecularFormula": "C6H12O6",
  "molecularMultiplier": 6,
  "molecularExact": true,
  "measuredMolarMass": 180.156,
  "impliedMolarMass": 180.156,
  "molarMassErrorPercent": 0
}
```

What it does deliberately:

- **The smallest ratio is not automatically 1** — `Fe` 111.69 g + `O` 47.997 g gives ratios 1 : 1.5, so
  subscripts come from the smallest multiplier `k ∈ [1, 8]` whose scaled ratios are all within ±0.05 of
  integers (here `k = 2` → `Fe2O3`), and `multiplier` is reported so the step is visible.
- **A molecular formula needs a whole multiple** — `molarMass` is only accepted as a molecular formula
  when it is within 2 % of an integer multiple of the empirical mass; otherwise `molecularExact: false`
  and a `note` say how far off it is.
- **No silent rounding** — when no multiplier fits, the closest fit is returned with `exact: false` and a
  note, because an analysis that does not support integer subscripts is a measurement problem, not a
  formula.
- **Percent data must add up** — `mode: "percent"` requires the parts to sum to 100 ± 1; a 90 % total is
  an error naming the total, not something to normalise away.
- **Hill order** — C, then H, then alphabetical (alphabetical without carbon), so the same analysis always
  prints the same string regardless of input order.
- **Errors instead of guesses** — unknown symbol, a duplicate element (`H` and `h`), a non-positive value,
  more than 12 elements, or a non-positive `molarMass` all return `valid: false` with the reason.

### `equilibrium` — weak acid / weak base ICE arithmetic

The calculation everyone skips: a weak acid does not dissociate completely, so `x` is the positive root
of `x² + K·x − K·C = 0`, not `√(K·C)`. This tool solves the quadratic exactly and also tells you what the
shortcut would have said.

```json
{ "mode": "weak_acid", "k": 1.8e-5, "concentration": 0.1 }
```

```json
{
  "valid": true,
  "mode": "weak_acid",
  "constantK": 1.8e-5,
  "pK": 4.7447274949,
  "concentration": 0.1,
  "x": 0.00133267097308,
  "hConc": 0.00133267097308,
  "ohConc": 7.50372762821e-12,
  "pH": 2.87527706155,
  "pOH": 11.1247229384,
  "degreeOfDissociation": 0.0133267097308,
  "percentDissociation": 1.33267097308,
  "approximation": { "x": 0.0013416407865, "errorPercent": 0.673070367938 },
  "kFromResult": 1.8e-5
}
```

- **The shortcut is reported, not assumed valid** — `approximation.errorPercent` is how much `√(K·C)`
  overstates `x`. It grows as the solution gets more dilute (0.21 % at 1 M acetic acid, 2.14 % at 0.01 M,
  6.93 % at 0.001 M ammonia with Kb = 1.8e-5); past 5 % a `note` says to use the quadratic answer.
- **Cancellation-free root** — `x = 2KC/(K + √(K²+4KC))` instead of `(−K + √(K²+4KC))/2`: algebraically the
  same root, but the second form loses precision exactly in the weak-acid regime (K ≪ C).
- **A back-check is part of the answer** — `kFromResult` recomputes `x²/(C − x)`, which must equal the K
  you supplied; if it does not, the numbers are wrong.
- **`weak_base` mirrors it** — Kb solves for `[OH-]`, then `pH` comes from Kw = 1e-14 at 25 °C.
- **Scope, stated plainly** — monoprotic acids/bases, no activity coefficients, no ionic strength, no
  polyprotic steps, no temperature dependence.

---
### `ksp_solubility` — solubility product → molar solubility

`√Ksp` is the easy part. The hard part is which ions the formula breaks into: `Ca(OH)2` is one Ca²⁺ and
two OH⁻ (so `Ksp = (s)(2s)² = 4s³`), and `Ca3(PO4)2` is `Ksp = (3s)³(2s)² = 108·s⁵` — a different
exponent *and* a different constant. Getting the grouping wrong changes the answer silently.

```json
{ "formula": "AgCl", "ksp": 1.8e-10 }
```

```json
{
  "valid": true,
  "formula": "AgCl",
  "ksp": 1.8e-10,
  "pKsp": 9.744727495,
  "cation": { "symbol": "Ag^+", "charge": 1, "count": 1, "concentration": 1.3416407864998738e-5 },
  "anion": { "symbol": "Cl^-", "charge": -1, "count": 1, "concentration": 1.3416407864998738e-5 },
  "stoichiometricFactor": 1,
  "ionCount": 2,
  "expression": "Ksp = [Ag^+]·[Cl^-] = s·s = s^2",
  "molarMass": 143.32,
  "molarSolubility": 1.3416407864998738e-5,
  "massSolubility": 0.00192283957521,
  "chargeBalance": 0,
  "notes": ["Ksp is used exactly as supplied; it is temperature dependent and this tool does not interpolate it"]
}
```

- **The grouping is verified, not guessed** — the formula is split into a cation half and an anion half,
  each half must match a known ion (composition *and* charge), and the candidate must be electroneutral.
  A formula that admits two readings is **refused with both candidates listed**, not resolved by
  preference: `Hg2Cl2` is `1 Hg2^2+ + 2 Cl^-` in solution chemistry and `2 Hg^+ + 2 Cl^-` in elementary
  texts, and those give different answers (`4s³` vs `16s⁴`). Pass `cation` (e.g. `"Hg2^2+"`) to choose.
- **The expression is printed** — `Ksp = [Ca^2+]^3·[PO4^3-]^2 = (3s)·(2s) = 108·s^5` is in the output so a
  reader can check the exponent and the constant without redoing the algebra.
- **Common-ion effect is solved, not approximated** — give `commonIon` + `commonIonConcentration` and the
  shifted equilibrium `(a·s + C)^a (b·s)^b = Ksp` is solved numerically (bisection), with the residual of
  its own back-substitution reported (`residual`, must be ~1e-16) and the reduction factor next to the
  plain-water solubility. AgCl in 0.01 M Cl⁻: 1.34e-5 → 1.80e-8 mol/L, 745× less.
- **Refusals name the reason** — hydrates (`Ksp` describes the anhydrous solid), charged or
  coefficient-prefixed formulas, mixed-valence oxides like `Fe3O4` that no ion pair reconstructs, and an
  ion that is not in the table all come back as `valid: false` with the reason and, where useful, the
  override hint.
- **Scope, stated plainly** — ideal solutions: no activity coefficients, no ionic strength, no hydrolysis,
  no complex formation, no temperature interpolation.

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
