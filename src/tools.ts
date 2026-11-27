/**
 * Tool definitions for dsh-chem: four deterministic chemistry tools exposed to
 * every agent via defineTool. Strict JSON-schema parameter surfaces, explicit
 * result interfaces matching the inferred output types.
 *
 * @module dsh-chem/tools
 */

import { defineTool, type ToolDefinition } from '@deepseek-ai/dsh-tools'
import { balanceEquation, renderBalancedEquation } from './balance.ts'
import { compositionToFormula } from './formula.ts'
import { weakEquilibrium } from './equil.ts'
import { kspSolubility } from './ksp.ts'
import { stoichiometry, type StoichResult } from './stoich.ts'
import {
  idealGas,
  gasDensity,
  molarMassFromDensity,
  molarVolume,
} from './gaslaw.ts'
import {
  AVOGADRO,
  bufferPh,
  formulaMass,
  parseFormula,
  phFromH,
  phFromOH,
  phFromPH,
  phFromPOH,
  resolveMolarMass,
  solutionPercent,
  type FormulaParse,
  type MassResult,
  type PhResult,
} from './core.ts'

export interface ToolSet {
  formula_mass: ToolDefinition
  mass_moles: ToolDefinition
  solution: ToolDefinition
  ph: ToolDefinition
  balance_equation: ToolDefinition
  gas_law: ToolDefinition
  stoichiometry: ToolDefinition
  formula_from_composition: ToolDefinition
  equilibrium: ToolDefinition
  ksp_solubility: ToolDefinition
}

/** Output-facing shape of `stoichiometry` (undefined keys are never emitted). */
export interface StoichiometryOut {
  valid: boolean
  equation?: string
  extent?: number
  limiting?: string[]
  reactants?: {
    formula: string
    coefficient: number
    molarMass: number
    availableMoles: number
    availableGrams: number
    consumedMoles: number
    leftoverMoles: number
    leftoverGrams: number
    limiting: boolean
  }[]
  products?: {
    formula: string
    coefficient: number
    molarMass: number
    moles: number
    grams: number
  }[]
  yieldProduct?: string
  theoreticalYieldMol?: number
  theoreticalYieldG?: number
  actualYieldG?: number
  percentYield?: number
  error?: string
}

export interface FormulaMassResult {
  valid: boolean
  formula?: string
  molarMass?: number
  elements?: {
    symbol: string
    name: string
    count: number
    atomicMass: number
    massContribution: number
    massFraction: number
  }[]
  hydrate?: { symbol: string; count: number }[]
  charge?: string
  coefficient?: number
  usedMassNumber?: boolean
  error?: string
}

export interface MassMolesResult {
  valid: boolean
  mode?: string
  moles?: number
  massG?: number
  particles?: number
  molarMass?: number
  formula?: string
  avogadro?: number
  error?: string
}

export interface SolutionResult {
  valid: boolean
  operation?: string
  molarity?: number
  massG?: number
  volumeMl?: number
  molarMass?: number
  formula?: string
  c1?: number
  v1?: number
  c2?: number
  v2?: number
  solved?: string
  percent?: number
  percentType?: string
  error?: string
}

export interface PhResultOut {
  valid: boolean
  operation?: string
  pH?: number
  pOH?: number
  hConc?: number
  ohConc?: number
  bufferRatio?: number
  pKa?: number
  error?: string
}

export interface BalanceTermOut {
  formula: string
  coefficient: number
}

export interface BalanceEquationResult {
  valid: boolean
  /** false when the species cannot balance at all (valid stays true: the input was understood). */
  balanced?: boolean
  /** false when several independent balances exist; the smallest one is returned. */
  unique?: boolean
  equation?: string
  reactants?: BalanceTermOut[]
  products?: BalanceTermOut[]
  elements?: { symbol: string; left: number; right: number }[]
  /** true when at least one species carried a charge (charge conservation checked too). */
  chargeBalanced?: boolean
  error?: string
}

export interface GasLawResult {
  valid: boolean
  operation?: string
  solved?: string
  pressureAtm?: number
  pressureKpa?: number
  pressureMmhg?: number
  volumeL?: number
  volumeMl?: number
  moles?: number
  temperatureK?: number
  temperatureC?: number
  molarVolumeLPerMol?: number
  molarVolumeMlPerMol?: number
  densityGPerL?: number
  densityGPerMl?: number
  molarMass?: number
  formula?: string
  r?: number
  error?: string
}

export interface CompositionOut {
  valid: boolean
  mode?: string
  empiricalFormula?: string
  empiricalMolarMass?: number
  multiplier?: number
  maxDeviation?: number
  exact?: boolean
  elements?: {
    symbol: string
    massG: number
    atomicMass: number
    moles: number
    ratio: number
    subscript: number
    massFraction: number
  }[]
  molecularFormula?: string
  molecularMultiplier?: number
  molecularExact?: boolean
  measuredMolarMass?: number
  impliedMolarMass?: number
  molarMassErrorPercent?: number
  note?: string
  error?: string
}

export interface EquilibriumOut {
  valid: boolean
  mode?: string
  constantK?: number
  pK?: number
  concentration?: number
  x?: number
  hConc?: number
  ohConc?: number
  pH?: number
  pOH?: number
  degreeOfDissociation?: number
  percentDissociation?: number
  approximation?: { x: number; errorPercent: number }
  kFromResult?: number
  note?: string
  error?: string
}

/** Output-facing shape of `ksp_solubility` (undefined keys are never emitted). */
export interface KspOut {
  valid: boolean
  formula?: string
  ksp?: number
  pKsp?: number
  cation?: { symbol: string; charge: number; count: number; concentration: number }
  anion?: { symbol: string; charge: number; count: number; concentration: number }
  stoichiometricFactor?: number
  ionCount?: number
  expression?: string
  molarMass?: number
  molarSolubility?: number
  massSolubility?: number
  chargeBalance?: number
  commonIon?: {
    symbol: string
    concentration: number
    molarSolubility: number
    massSolubility: number
    reductionFactor: number
    residual: number
  }
  notes?: string[]
  candidates?: string[]
  error?: string
}

function fmt(value: number): string {
  if (!Number.isFinite(value)) return String(value)
  if (value !== 0 && (Math.abs(value) >= 1e7 || Math.abs(value) < 1e-4)) {
    return value.toExponential(6)
  }
  return String(Number(value.toPrecision(8)))
}

/* ------------------------------------------------------------------ */
/* render helpers                                                      */
/* ------------------------------------------------------------------ */

function renderFormulaMass(value: unknown): string {
  const v = value as FormulaMassResult
  if (!v.valid) return `formula_mass failed: ${v.error}`
  const lines: string[] = [`${v.formula}: molar mass ${fmt(v.molarMass!)} g/mol`]
  for (const e of v.elements!) {
    lines.push(`  ${e.symbol.padEnd(3)} × ${String(e.count).padStart(3)} = ${fmt(e.massContribution)} g/mol (${fmt(e.massFraction * 100)} %)`)
  }
  if (v.hydrate && v.hydrate.length > 0) {
    const water = v.hydrate.map((h) => `${h.symbol}${h.count}`).join('')
    lines.push(`  hydrate water: ${water}`)
  }
  if (v.charge) lines.push(`  ionic charge: ${v.charge}`)
  if (v.usedMassNumber) lines.push('  note: some elements use most-stable-isotope mass numbers (radioelements)')
  return lines.join('\n')
}

function renderMassMoles(value: unknown): string {
  const v = value as MassMolesResult
  if (!v.valid) return `mass_moles failed: ${v.error}`
  const lines: string[] = [`mode: ${v.mode}`]
  if (v.molarMass !== undefined) lines.push(`molar mass: ${fmt(v.molarMass)} g/mol`)
  if (v.massG !== undefined) lines.push(`mass: ${fmt(v.massG)} g`)
  if (v.moles !== undefined) lines.push(`amount: ${fmt(v.moles)} mol`)
  if (v.particles !== undefined) lines.push(`particles: ${fmt(v.particles)}`)
  return lines.join('\n')
}

function renderSolution(value: unknown): string {
  const v = value as SolutionResult
  if (!v.valid) return `solution failed: ${v.error}`
  const lines: string[] = [`operation: ${v.operation}`]
  if (v.molarity !== undefined) lines.push(`molarity: ${fmt(v.molarity)} mol/L`)
  if (v.massG !== undefined) lines.push(`solute mass: ${fmt(v.massG)} g`)
  if (v.volumeMl !== undefined) lines.push(`volume: ${fmt(v.volumeMl)} mL`)
  if (v.percent !== undefined) lines.push(`${v.percentType} percent: ${fmt(v.percent)} %`)
  if (v.solved) lines.push(`solved: ${v.solved}`)
  if (v.c1 !== undefined && v.c2 !== undefined) lines.push(`dilution: ${fmt(v.c1)} M × ${fmt(v.v1!)} mL = ${fmt(v.c2)} M × ${fmt(v.v2!)} mL`)
  return lines.join('\n')
}

function renderBalance(value: unknown): string {
  const v = value as BalanceEquationResult
  if (!v.valid) return `balance_equation failed: ${v.error}`
  if (!v.balanced) return `balance_equation: cannot balance — ${v.error}`
  const lines: string[] = [v.equation!]
  lines.push(v.unique === false
    ? 'note: several independent balances exist — this is the smallest one by total coefficient sum'
    : 'balanced (unique smallest integer solution)')
  for (const e of v.elements!) {
    lines.push(`  ${e.symbol.padEnd(3)} left ${String(e.left).padStart(3)} = right ${String(e.right).padStart(3)}`)
  }
  if (v.chargeBalanced) lines.push('  ionic charge is conserved as well')
  return lines.join('\n')
}

function renderGasLaw(value: unknown): string {
  const v = value as GasLawResult
  if (!v.valid) return `gas_law failed: ${v.error}`
  const lines: string[] = [`operation: ${v.operation}`]
  if (v.solved) lines.push(`solved for: ${v.solved}`)
  if (v.formula) lines.push(`formula: ${v.formula}`)
  if (v.pressureAtm !== undefined) lines.push(`pressure: ${fmt(v.pressureAtm)} atm (${fmt(v.pressureKpa!)} kPa, ${fmt(v.pressureMmhg!)} mmHg)`)
  if (v.volumeL !== undefined) lines.push(`volume: ${fmt(v.volumeL)} L (${fmt(v.volumeMl!)} mL)`)
  if (v.moles !== undefined) lines.push(`amount: ${fmt(v.moles)} mol`)
  if (v.temperatureK !== undefined) lines.push(`temperature: ${fmt(v.temperatureK)} K (${fmt(v.temperatureC!)} °C)`)
  if (v.molarVolumeLPerMol !== undefined) lines.push(`molar volume: ${fmt(v.molarVolumeLPerMol)} L/mol (${fmt(v.molarVolumeMlPerMol!)} mL/mol)`)
  if (v.densityGPerL !== undefined) lines.push(`density: ${fmt(v.densityGPerL)} g/L`)
  if (v.molarMass !== undefined) lines.push(`molar mass: ${fmt(v.molarMass)} g/mol`)
  lines.push(`R = ${fmt(v.r!)} L·atm/(mol·K)`)
  return lines.join('\n')
}

function renderStoichiometry(value: unknown): string {
  const v = value as StoichiometryOut
  if (!v.valid) return `stoichiometry failed: ${v.error}`
  const lines: string[] = [v.equation!]
  lines.push(`extent of reaction: ${fmt(v.extent!)} mol`)
  lines.push(`limiting reagent: ${v.limiting!.join(' + ')}`)
  for (const r of v.reactants!) {
    lines.push(`  ${r.formula.padEnd(10)} × ${String(r.coefficient).padStart(2)}  available ${fmt(r.availableMoles)} mol (${fmt(r.availableGrams)} g)`
      + `  consumed ${fmt(r.consumedMoles)} mol  left ${fmt(r.leftoverMoles)} mol (${fmt(r.leftoverGrams)} g)`
      + (r.limiting ? '  [limiting]' : ''))
  }
  for (const p of v.products!) {
    lines.push(`product ${p.formula} × ${p.coefficient}: ${fmt(p.moles)} mol (${fmt(p.grams)} g)`)
  }
  if (v.yieldProduct !== undefined) {
    lines.push(`theoretical yield of ${v.yieldProduct}: ${fmt(v.theoreticalYieldG!)} g (${fmt(v.theoreticalYieldMol!)} mol)`)
  }
  if (v.percentYield !== undefined) {
    lines.push(`percent yield: ${fmt(v.percentYield)} % (measured ${fmt(v.actualYieldG!)} g)`)
  }
  return lines.join('\n')
}

function renderComposition(value: unknown): string {
  const v = value as CompositionOut
  if (!v.valid) return `formula_from_composition failed: ${v.error}`
  const lines: string[] = [`empirical formula: ${v.empiricalFormula}  (M = ${fmt(v.empiricalMolarMass!)} g/mol)`]
  for (const e of v.elements!) {
    lines.push(`  ${e.symbol.padEnd(3)} ${fmt(e.massG).padStart(10)} g  ${fmt(e.moles).padStart(10)} mol`
      + `  ratio ${fmt(e.ratio).padStart(9)}  → ${e.subscript}  (${fmt(e.massFraction)} % of the formula)`)
  }
  lines.push(`integer multiplier: ${v.multiplier} (max deviation ${fmt(v.maxDeviation!)})`)
  if (v.exact === false) lines.push('⚠️ the analysis does not support integer subscripts within tolerance')
  if (v.molecularFormula !== undefined) {
    lines.push(`molecular formula: ${v.molecularFormula}  (×${v.molecularMultiplier} of the empirical formula,`
      + ` implied M = ${fmt(v.impliedMolarMass!)} g/mol vs measured ${fmt(v.measuredMolarMass!)} g/mol,`
      + ` ${fmt(v.molarMassErrorPercent!)} % off)`)
    if (v.molecularExact === false) lines.push('⚠️ the measured molar mass is not a whole multiple of the empirical mass')
  }
  if (v.note) lines.push(v.note)
  return lines.join('\n')
}

function renderEquilibrium(value: unknown): string {
  const v = value as EquilibriumOut
  if (!v.valid) return `equilibrium failed: ${v.error}`
  const label = v.mode === 'weak_base' ? 'Kb' : 'Ka'
  const lines: string[] = [`${v.mode}: ${label} = ${fmt(v.constantK!)} (p${label} = ${fmt(v.pK!)}), C = ${fmt(v.concentration!)} mol/L`]
  lines.push(`dissociated: x = ${fmt(v.x!)} mol/L (${fmt(v.percentDissociation!)} % of C)`)
  lines.push(`[H+] = ${fmt(v.hConc!)} mol/L   [OH-] = ${fmt(v.ohConc!)} mol/L`)
  lines.push(`pH = ${fmt(v.pH!)}   pOH = ${fmt(v.pOH!)}`)
  lines.push(`√(K·C) shortcut would give x = ${fmt(v.approximation!.x)} (${fmt(v.approximation!.errorPercent)} % high)`)
  lines.push(`back-check: x²/(C − x) = ${fmt(v.kFromResult!)} (must equal ${label})`)
  if (v.note) lines.push(v.note)
  return lines.join('\n')
}

function renderKsp(value: unknown): string {
  const v = value as KspOut
  if (!v.valid) {
    const lines = [`ksp_solubility failed: ${v.error}`]
    if (v.candidates && v.candidates.length > 0) lines.push(`candidates: ${v.candidates.join(' | ')}`)
    return lines.join('\n')
  }
  const lines: string[] = [
    `${v.formula}: ${v.cation!.count} ${v.cation!.symbol} + ${v.anion!.count} ${v.anion!.symbol}`
    + ` (charge balance ${fmt(v.chargeBalance!)}), Ksp = ${fmt(v.ksp!)} (pKsp = ${fmt(v.pKsp!)})`,
  ]
  lines.push(v.expression!)
  lines.push(`molar solubility: ${fmt(v.molarSolubility!)} mol/L (${fmt(v.massSolubility!)} g/L, M = ${fmt(v.molarMass!)} g/mol)`)
  lines.push(`at saturation: [${v.cation!.symbol}] = ${fmt(v.cation!.concentration!)} mol/L, [${v.anion!.symbol}] = ${fmt(v.anion!.concentration!)} mol/L`)
  if (v.commonIon) {
    const c = v.commonIon
    lines.push(
      `common ion ${c.symbol} at ${fmt(c.concentration)} mol/L: solubility drops to ${fmt(c.molarSolubility)} mol/L `
      + `(${fmt(c.reductionFactor)}× less), ${fmt(c.massSolubility)} g/L`,
    )
  }
  for (const note of v.notes ?? []) lines.push(note)
  return lines.join('\n')
}

function renderPh(value: unknown): string {
  const v = value as PhResultOut
  if (!v.valid) return `ph failed: ${v.error}`
  const lines: string[] = [`operation: ${v.operation}`]
  if (v.pH !== undefined) lines.push(`pH: ${fmt(v.pH)}`)
  if (v.pOH !== undefined) lines.push(`pOH: ${fmt(v.pOH)}`)
  if (v.hConc !== undefined) lines.push(`[H+]: ${fmt(v.hConc)} mol/L`)
  if (v.ohConc !== undefined) lines.push(`[OH-]: ${fmt(v.ohConc)} mol/L`)
  if (v.bufferRatio !== undefined) lines.push(`[A-]/[HA] ratio: ${fmt(v.bufferRatio)}`)
  return lines.join('\n')
}

/* ------------------------------------------------------------------ */
/* tool definitions                                                    */
/* ------------------------------------------------------------------ */

export function buildChemTools(): ToolSet {
  const formula_mass = defineTool({
    name: 'formula_mass',
    description: 'Parse a chemical formula (nested parentheses, hydrates like "CuSO4·5H2O", ionic charges '
      + 'like "Fe3+" or "SO4^2-", leading coefficients like "2H2O") and compute its exact molar mass with a '
      + 'per-element breakdown and mass percentages. Uses IUPAC standard atomic weights; radioelements use '
      + 'most-stable-isotope mass numbers. Use this instead of mental arithmetic for molar-mass math — '
      + 'symbols are case-sensitive ("Co" is cobalt, "CO" is carbon monoxide).',
    parameters: {
      formula: { type: 'string', required: true, description: 'Chemical formula, e.g. "H2SO4", "Ca(OH)2", "CuSO4·5H2O", "Fe2(SO4)3", "CH3CH2OH". Accepts "·" or "*" as hydrate separator.' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          formula: { type: 'string' },
          molarMass: { type: 'number' },
          elements: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                symbol: { type: 'string', required: true },
                name: { type: 'string', required: true },
                count: { type: 'number', required: true },
                atomicMass: { type: 'number', required: true },
                massContribution: { type: 'number', required: true },
                massFraction: { type: 'number', required: true },
              },
            },
          },
          hydrate: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                symbol: { type: 'string', required: true },
                count: { type: 'number', required: true },
              },
            },
          },
          charge: { type: 'string' },
          coefficient: { type: 'number' },
          usedMassNumber: { type: 'boolean' },
          error: { type: 'string' },
        },
      },
      render: (_args: { formula: string }, value: unknown) => [{ type: 'text', text: renderFormulaMass(value) }],
    },
    async execute(args: { formula: string }): Promise<FormulaMassResult> {
      const parsed: FormulaParse = parseFormula(args.formula)
      if (!parsed.ok) return { valid: false, error: parsed.error }
      const mass: MassResult = formulaMass(parsed)
      const result: FormulaMassResult = {
        valid: true,
        formula: parsed.normalized,
        molarMass: mass.molarMass,
        elements: mass.breakdown.map((b) => ({
          symbol: b.symbol,
          name: b.name,
          count: b.count,
          atomicMass: b.atomicMass,
          massContribution: b.massContribution,
          massFraction: b.massFraction,
        })),
        usedMassNumber: mass.usedMassNumber,
      }
      if (parsed.hydrate && parsed.hydrate.length > 0) result.hydrate = parsed.hydrate
      if (parsed.charge) result.charge = parsed.charge
      if (parsed.coefficient !== null) result.coefficient = parsed.coefficient
      return result
    },
  })

  const mass_moles = defineTool({
    name: 'mass_moles',
    description: 'Convert between mass (g), amount of substance (mol) and particle count using the molar mass '
      + '(from a formula or given directly) and the Avogadro constant 6.02214076e23. Six modes: mass_to_moles, '
      + 'moles_to_mass, mass_to_particles, particles_to_mass, moles_to_particles, particles_to_moles. '
      + 'Use this instead of mental arithmetic for mole conversions.',
    parameters: {
      mode: {
        type: 'string',
        required: true,
        enum: ['mass_to_moles', 'moles_to_mass', 'mass_to_particles', 'particles_to_mass', 'moles_to_particles', 'particles_to_moles'],
        description: 'Conversion direction.',
      },
      value: { type: 'number', required: true, description: 'The quantity to convert: grams for mass_* modes, moles for moles_* modes, particle count for particles_* modes.' },
      formula: { type: 'string', description: 'Chemical formula whose molar mass to use, e.g. "H2SO4". Takes precedence over molarMass when both are given.' },
      molarMass: { type: 'number', description: 'Molar mass in g/mol (used when formula is not given).' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          mode: { type: 'string' },
          moles: { type: 'number' },
          massG: { type: 'number' },
          particles: { type: 'number' },
          molarMass: { type: 'number' },
          formula: { type: 'string' },
          avogadro: { type: 'number' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderMassMoles(value) }],
    },
    async execute(args: { mode: string; value: number; formula?: string; molarMass?: number }): Promise<MassMolesResult> {
      const mode = args.mode
      if (!Number.isFinite(args.value) || args.value <= 0) {
        return { valid: false, error: `value must be a positive number (got ${args.value})` }
      }
      const needsMass = mode === 'mass_to_moles' || mode === 'mass_to_particles' || mode === 'particles_to_mass' || mode === 'moles_to_mass'
      if (needsMass) {
        const resolved = resolveMolarMass(args.formula, args.molarMass)
        if (resolved.error) return { valid: false, error: resolved.error }
        const M = resolved.mass
        const result: MassMolesResult = { valid: true, mode, molarMass: M, avogadro: AVOGADRO }
        if (resolved.formula) result.formula = resolved.formula
        if (mode === 'mass_to_moles') result.moles = args.value / M
        else if (mode === 'mass_to_particles') result.particles = (args.value / M) * AVOGADRO
        else if (mode === 'particles_to_mass') result.massG = (args.value / AVOGADRO) * M
        else result.massG = args.value * M
        return result
      }
      const result: MassMolesResult = { valid: true, mode, avogadro: AVOGADRO }
      if (mode === 'moles_to_particles') result.particles = args.value * AVOGADRO
      else if (mode === 'particles_to_moles') result.moles = args.value / AVOGADRO
      else return { valid: false, error: `unknown mode "${mode}" (expected one of mass_to_moles, moles_to_mass, mass_to_particles, particles_to_mass, moles_to_particles, particles_to_moles)` }
      return result
    },
  })

  const solution = defineTool({
    name: 'solution',
    description: 'Solution concentration math. Four operations: "molarity" (mass of solute + volume -> mol/L), '
      + '"prepare" (target molarity + volume -> solute mass to weigh), "dilute" (C1V1 = C2V2, give any three of '
      + 'c1/v1/c2/v2 in M and mL and the missing one is computed), and "percent" (w/v, v/v or mass percent). '
      + 'Use this instead of mental arithmetic for dilution and concentration problems.',
    parameters: {
      operation: { type: 'string', required: true, enum: ['molarity', 'prepare', 'dilute', 'percent'], description: 'Which calculation to perform.' },
      massG: { type: 'number', description: 'Solute mass in grams (molarity operation).' },
      volumeMl: { type: 'number', description: 'Solution volume in milliliters (molarity / prepare operations).' },
      molarity: { type: 'number', description: 'Molarity in mol/L (prepare operation).' },
      formula: { type: 'string', description: 'Solute chemical formula, e.g. "NaCl". Takes precedence over molarMass.' },
      molarMass: { type: 'number', description: 'Solute molar mass in g/mol (used when formula is not given).' },
      c1: { type: 'number', description: 'Stock concentration in mol/L (dilute operation).' },
      v1: { type: 'number', description: 'Stock volume in mL (dilute operation).' },
      c2: { type: 'number', description: 'Final concentration in mol/L (dilute operation).' },
      v2: { type: 'number', description: 'Final volume in mL (dilute operation).' },
      percentType: { type: 'string', enum: ['wv', 'vv', 'mass'], description: 'Percent kind: w/v (g per 100 mL), v/v (mL per 100 mL), mass (mass %).' },
      soluteG: { type: 'number', description: 'Solute grams (percent w/v).' },
      solutionMl: { type: 'number', description: 'Solution milliliters (percent w/v).' },
      soluteMl: { type: 'number', description: 'Solute milliliters (percent v/v).' },
      totalMl: { type: 'number', description: 'Total solution milliliters (percent v/v).' },
      soluteMassG: { type: 'number', description: 'Solute mass in grams (percent mass).' },
      totalMassG: { type: 'number', description: 'Total mixture mass in grams (percent mass).' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          operation: { type: 'string' },
          molarity: { type: 'number' },
          massG: { type: 'number' },
          volumeMl: { type: 'number' },
          molarMass: { type: 'number' },
          formula: { type: 'string' },
          c1: { type: 'number' },
          v1: { type: 'number' },
          c2: { type: 'number' },
          v2: { type: 'number' },
          solved: { type: 'string' },
          percent: { type: 'number' },
          percentType: { type: 'string' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderSolution(value) }],
    },
    async execute(args: {
      operation: string
      massG?: number
      volumeMl?: number
      molarity?: number
      formula?: string
      molarMass?: number
      c1?: number
      v1?: number
      c2?: number
      v2?: number
      percentType?: string
      soluteG?: number
      solutionMl?: number
      soluteMl?: number
      totalMl?: number
      soluteMassG?: number
      totalMassG?: number
    }): Promise<SolutionResult> {
      const op = args.operation
      if (op === 'molarity' || op === 'prepare') {
        const resolved = resolveMolarMass(args.formula, args.molarMass)
        if (resolved.error) return { valid: false, error: resolved.error }
        if (op === 'molarity') {
          if (args.massG === undefined || args.volumeMl === undefined) {
            return { valid: false, error: 'operation "molarity" requires massG (grams) and volumeMl (milliliters)' }
          }
          if (!Number.isFinite(args.massG) || args.massG <= 0) return { valid: false, error: `massG must be positive (got ${args.massG})` }
          if (!Number.isFinite(args.volumeMl) || args.volumeMl <= 0) return { valid: false, error: `volumeMl must be positive (got ${args.volumeMl})` }
          const result: SolutionResult = {
            valid: true,
            operation: 'molarity',
            molarity: args.massG / resolved.mass / (args.volumeMl / 1000),
            massG: args.massG,
            volumeMl: args.volumeMl,
            molarMass: resolved.mass,
          }
          if (resolved.formula) result.formula = resolved.formula
          return result
        }
        if (args.molarity === undefined || args.volumeMl === undefined) {
          return { valid: false, error: 'operation "prepare" requires molarity (mol/L) and volumeMl (milliliters)' }
        }
        if (!Number.isFinite(args.molarity) || args.molarity <= 0) return { valid: false, error: `molarity must be positive (got ${args.molarity})` }
        if (!Number.isFinite(args.volumeMl) || args.volumeMl <= 0) return { valid: false, error: `volumeMl must be positive (got ${args.volumeMl})` }
        const result: SolutionResult = {
          valid: true,
          operation: 'prepare',
          massG: args.molarity * (args.volumeMl / 1000) * resolved.mass,
          molarity: args.molarity,
          volumeMl: args.volumeMl,
          molarMass: resolved.mass,
        }
        if (resolved.formula) result.formula = resolved.formula
        return result
      }
      if (op === 'dilute') {
        const given = [args.c1 !== undefined, args.v1 !== undefined, args.c2 !== undefined, args.v2 !== undefined]
        const count = given.filter(Boolean).length
        if (count !== 3) {
          return { valid: false, error: 'operation "dilute" requires exactly three of c1, v1, c2, v2 (C1V1 = C2V2); give the missing one as the unknown' }
        }
        const values: Record<'c1' | 'v1' | 'c2' | 'v2', number | undefined> = { c1: args.c1, v1: args.v1, c2: args.c2, v2: args.v2 }
        for (const key of ['c1', 'v1', 'c2', 'v2'] as const) {
          const val = values[key]
          if (val !== undefined && (!Number.isFinite(val) || val <= 0)) {
            return { valid: false, error: `${key} must be positive (got ${val})` }
          }
        }
        let solvedKey: 'c1' | 'v1' | 'c2' | 'v2'
        let solvedValue: number
        if (args.c1 === undefined) {
          solvedKey = 'c1'
          solvedValue = (args.c2! * args.v2!) / args.v1!
        } else if (args.v1 === undefined) {
          solvedKey = 'v1'
          solvedValue = (args.c2! * args.v2!) / args.c1
        } else if (args.c2 === undefined) {
          solvedKey = 'c2'
          solvedValue = (args.c1 * args.v1!) / args.v2!
        } else {
          solvedKey = 'v2'
          solvedValue = (args.c1 * args.v1!) / args.c2!
        }
        const result: SolutionResult = {
          valid: true,
          operation: 'dilute',
          c1: args.c1 === undefined ? solvedValue : args.c1,
          v1: args.v1 === undefined ? solvedValue : args.v1,
          c2: args.c2 === undefined ? solvedValue : args.c2,
          v2: args.v2 === undefined ? solvedValue : args.v2,
          solved: solvedKey,
        }
        return result
      }
      if (op === 'percent') {
        const pct = solutionPercent({
          percentType: (args.percentType ?? 'wv') as 'wv' | 'vv' | 'mass',
          soluteG: args.soluteG,
          solutionMl: args.solutionMl,
          soluteMl: args.soluteMl,
          totalMl: args.totalMl,
          soluteMassG: args.soluteMassG,
          totalMassG: args.totalMassG,
        })
        if (!pct.ok) return { valid: false, error: pct.error }
        return { valid: true, operation: 'percent', percent: pct.percent, percentType: args.percentType ?? 'wv' }
      }
      return { valid: false, error: `unknown operation "${op}" (expected molarity, prepare, dilute or percent)` }
    },
  })

  const ph = defineTool({
    name: 'ph',
    description: 'pH chemistry at 25 °C (Kw = 1e-14). Five operations: "from_h" ([H+] in mol/L -> pH/pOH), '
      + '"from_oh" ([OH-] -> pH/pOH), "from_ph" (pH -> [H+]), "from_poh" (pOH -> [OH-]), and "buffer" '
      + '(Henderson-Hasselbalch: pH = pKa + log10([A-]/[HA]) from pKa, base and acid concentrations). '
      + 'Use this instead of mental log/antilog arithmetic for pH problems.',
    parameters: {
      operation: { type: 'string', required: true, enum: ['from_h', 'from_oh', 'from_ph', 'from_poh', 'buffer'], description: 'Which pH calculation to perform.' },
      hConc: { type: 'number', description: '[H+] in mol/L (from_h operation), must be > 0.' },
      ohConc: { type: 'number', description: '[OH-] in mol/L (from_oh operation), must be > 0.' },
      pH: { type: 'number', description: 'pH value (from_ph operation), typically 0-14 but any value in [-2, 16] is accepted.' },
      pOH: { type: 'number', description: 'pOH value (from_poh operation).' },
      pKa: { type: 'number', description: 'Acid dissociation constant as pKa (buffer operation), e.g. 4.76 for acetic acid.' },
      acidConc: { type: 'number', description: '[HA] weak-acid concentration in mol/L (buffer operation), must be > 0.' },
      baseConc: { type: 'number', description: '[A-] conjugate-base concentration in mol/L (buffer operation), must be > 0.' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          operation: { type: 'string' },
          pH: { type: 'number' },
          pOH: { type: 'number' },
          hConc: { type: 'number' },
          ohConc: { type: 'number' },
          bufferRatio: { type: 'number' },
          pKa: { type: 'number' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderPh(value) }],
    },
    async execute(args: {
      operation: string
      hConc?: number
      ohConc?: number
      pH?: number
      pOH?: number
      pKa?: number
      acidConc?: number
      baseConc?: number
    }): Promise<PhResultOut> {
      const op = args.operation
      const toResult = (r: PhResult): PhResultOut => ({
        valid: true,
        operation: op,
        pH: r.pH,
        pOH: r.pOH,
        hConc: r.hConc,
        ohConc: r.ohConc,
      })
      if (op === 'from_h') {
        if (args.hConc === undefined) return { valid: false, error: 'operation "from_h" requires hConc (mol/L)' }
        if (!Number.isFinite(args.hConc) || args.hConc <= 0) return { valid: false, error: `hConc must be positive (got ${args.hConc})` }
        return toResult(phFromH(args.hConc))
      }
      if (op === 'from_oh') {
        if (args.ohConc === undefined) return { valid: false, error: 'operation "from_oh" requires ohConc (mol/L)' }
        if (!Number.isFinite(args.ohConc) || args.ohConc <= 0) return { valid: false, error: `ohConc must be positive (got ${args.ohConc})` }
        return toResult(phFromOH(args.ohConc))
      }
      if (op === 'from_ph') {
        if (args.pH === undefined) return { valid: false, error: 'operation "from_ph" requires pH' }
        if (!Number.isFinite(args.pH) || args.pH < -2 || args.pH > 16) return { valid: false, error: `pH must be in [-2, 16] (got ${args.pH})` }
        return toResult(phFromPH(args.pH))
      }
      if (op === 'from_poh') {
        if (args.pOH === undefined) return { valid: false, error: 'operation "from_poh" requires pOH' }
        if (!Number.isFinite(args.pOH) || args.pOH < -2 || args.pOH > 16) return { valid: false, error: `pOH must be in [-2, 16] (got ${args.pOH})` }
        return toResult(phFromPOH(args.pOH))
      }
      if (op === 'buffer') {
        if (args.pKa === undefined || args.acidConc === undefined || args.baseConc === undefined) {
          return { valid: false, error: 'operation "buffer" requires pKa, acidConc and baseConc (all in mol/L except pKa)' }
        }
        if (!Number.isFinite(args.pKa)) return { valid: false, error: `pKa must be a number (got ${args.pKa})` }
        if (!Number.isFinite(args.acidConc) || args.acidConc <= 0) return { valid: false, error: `acidConc must be positive (got ${args.acidConc})` }
        if (!Number.isFinite(args.baseConc) || args.baseConc <= 0) return { valid: false, error: `baseConc must be positive (got ${args.baseConc})` }
        const buffered = bufferPh(args.pKa, args.baseConc, args.acidConc)
        return { valid: true, operation: 'buffer', pH: buffered.pH, bufferRatio: buffered.ratio, pKa: args.pKa }
      }
      return { valid: false, error: `unknown operation "${op}" (expected from_h, from_oh, from_ph, from_poh or buffer)` }
    },
  })

  const balance_equation = defineTool({
    name: 'balance_equation',
    description: 'Balance a chemical equation and return the smallest integer coefficients '
      + '(exact rational arithmetic — never floating point, so no 2.9999997 coefficients). '
      + 'Write species separated by " + " and one arrow ("->", "=" or "→"), e.g. '
      + '"C3H8 + O2 -> CO2 + H2O" or "MnO4- + Fe2+ + H+ -> Mn2+ + Fe3+ + H2O". '
      + 'Handles nested parentheses, hydrates (CuSO4·5H2O) and ionic charges; when any species '
      + 'carries a charge, charge conservation is enforced as an extra constraint. Per-element '
      + 'atom counts on both sides are returned for verification. Species may be written in the '
      + 'compact form ("Fe+O2->Fe2O3") but charged species must use the spaced form ("Na+ + Cl-") '
      + 'because a bare "+" next to a symbol is a charge, not a separator. Use this instead of '
      + 'eyeballing coefficients — it reports "cannot balance" rather than inventing a fake result.',
    parameters: {
      equation: { type: 'string', required: true, description: 'The unbalanced equation, e.g. "Fe + O2 -> Fe2O3", "C6H12O6 + O2 -> CO2 + H2O".' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          balanced: { type: 'boolean' },
          unique: { type: 'boolean' },
          equation: { type: 'string' },
          reactants: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                formula: { type: 'string', required: true },
                coefficient: { type: 'number', required: true },
              },
            },
          },
          products: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                formula: { type: 'string', required: true },
                coefficient: { type: 'number', required: true },
              },
            },
          },
          elements: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                symbol: { type: 'string', required: true },
                left: { type: 'number', required: true },
                right: { type: 'number', required: true },
              },
            },
          },
          chargeBalanced: { type: 'boolean' },
          error: { type: 'string' },
        },
      },
      render: (_args: { equation: string }, value: unknown) => [{ type: 'text', text: renderBalance(value) }],
    },
    async execute(args: { equation: string }): Promise<BalanceEquationResult> {
      const result = balanceEquation(args.equation)
      if (!result.ok) {
        if (result.kind === 'parse') return { valid: false, error: result.error }
        return { valid: true, balanced: false, error: result.error }
      }
      return {
        valid: true,
        balanced: true,
        unique: result.unique,
        equation: renderBalancedEquation(result.reactants, result.products),
        reactants: result.reactants,
        products: result.products,
        elements: result.elements,
        chargeBalanced: result.chargeBalanced,
      }
    },
  })

  const gas_law = defineTool({
    name: 'gas_law',
    description: 'Ideal-gas arithmetic with fixed units (pressure atm — also reported in kPa and mmHg; volume L — '
      + 'also mL; amount mol; temperature K — also °C). Four operations: "ideal" solves PV = nRT for the one '
      + 'quantity you leave out (give exactly three of pressureAtm, volumeL, moles, temperatureK); "molar_volume" '
      + 'gives V/n = RT/P (defaults 0 °C, 1 atm → 22.414 L/mol); "density" gives ρ = P·M/(R·T) from a formula or a '
      + 'given molarMass; "molar_mass" back-calculates M = ρ·R·T/P from a measured density. R = 8.31446261815324 / '
      + '101.325 L·atm/(mol·K). Use this instead of mental PV=nRT arithmetic and unit conversions.',
    parameters: {
      operation: { type: 'string', required: true, enum: ['ideal', 'molar_volume', 'density', 'molar_mass'], description: 'Which gas calculation to perform.' },
      pressureAtm: { type: 'number', description: 'Pressure in atm (must be > 0).' },
      volumeL: { type: 'number', description: 'Volume in litres (must be > 0).' },
      moles: { type: 'number', description: 'Amount of gas in mol (must be > 0).' },
      temperatureK: { type: 'number', description: 'Absolute temperature in K (must be > 0); omit to solve for it in "ideal".' },
      formula: { type: 'string', description: 'Gas formula for "density", e.g. "CO2". Takes precedence over molarMass.' },
      molarMass: { type: 'number', description: 'Molar mass in g/mol for "density" when no formula is given.' },
      densityGPerL: { type: 'number', description: 'Measured gas density in g/L for the "molar_mass" operation.' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          operation: { type: 'string' },
          solved: { type: 'string' },
          pressureAtm: { type: 'number' },
          pressureKpa: { type: 'number' },
          pressureMmhg: { type: 'number' },
          volumeL: { type: 'number' },
          volumeMl: { type: 'number' },
          moles: { type: 'number' },
          temperatureK: { type: 'number' },
          temperatureC: { type: 'number' },
          molarVolumeLPerMol: { type: 'number' },
          molarVolumeMlPerMol: { type: 'number' },
          densityGPerL: { type: 'number' },
          densityGPerMl: { type: 'number' },
          molarMass: { type: 'number' },
          formula: { type: 'string' },
          r: { type: 'number' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderGasLaw(value) }],
    },
    async execute(args: {
      operation: string
      pressureAtm?: number
      volumeL?: number
      moles?: number
      temperatureK?: number
      formula?: string
      molarMass?: number
      densityGPerL?: number
    }): Promise<GasLawResult> {
      const op = args.operation
      if (op === 'ideal') {
        const solved = idealGas({
          pressureAtm: args.pressureAtm,
          volumeL: args.volumeL,
          moles: args.moles,
          temperatureK: args.temperatureK,
        })
        if (!solved.ok) return { valid: false, operation: op, error: solved.error }
        return {
          valid: true,
          operation: op,
          solved: solved.solved,
          pressureAtm: solved.pressureAtm,
          pressureKpa: solved.pressureKpa,
          pressureMmhg: solved.pressureMmhg,
          volumeL: solved.volumeL,
          volumeMl: solved.volumeMl,
          moles: solved.moles,
          temperatureK: solved.temperatureK,
          temperatureC: solved.temperatureC,
          r: solved.r,
        }
      }
      if (op === 'molar_volume') {
        const solved = molarVolume(args.temperatureK ?? 273.15, args.pressureAtm ?? 1)
        if (!solved.ok) return { valid: false, operation: op, error: solved.error }
        return {
          valid: true,
          operation: op,
          pressureAtm: solved.pressureAtm,
          pressureKpa: solved.pressureKpa,
          temperatureK: solved.temperatureK,
          temperatureC: solved.temperatureC,
          molarVolumeLPerMol: solved.molarVolumeLPerMol,
          molarVolumeMlPerMol: solved.molarVolumeMlPerMol,
          r: solved.r,
        }
      }
      if (op === 'density') {
        const mass = resolveMolarMass(args.formula, args.molarMass)
        if (mass.error !== null) return { valid: false, operation: op, error: mass.error }
        if (args.pressureAtm === undefined) return { valid: false, operation: op, error: 'operation "density" requires pressureAtm' }
        if (args.temperatureK === undefined) return { valid: false, operation: op, error: 'operation "density" requires temperatureK' }
        const solved = gasDensity(args.pressureAtm, mass.mass, args.temperatureK)
        if (!solved.ok) return { valid: false, operation: op, error: solved.error }
        const out: GasLawResult = {
          valid: true,
          operation: op,
          pressureAtm: solved.pressureAtm,
          temperatureK: solved.temperatureK,
          temperatureC: solved.temperatureC,
          molarMass: solved.molarMass,
          densityGPerL: solved.densityGPerL,
          densityGPerMl: solved.densityGPerMl,
          r: solved.r,
        }
        if (mass.formula !== null) out.formula = mass.formula
        return out
      }
      if (op === 'molar_mass') {
        if (args.densityGPerL === undefined) return { valid: false, operation: op, error: 'operation "molar_mass" requires densityGPerL' }
        if (args.pressureAtm === undefined) return { valid: false, operation: op, error: 'operation "molar_mass" requires pressureAtm' }
        if (args.temperatureK === undefined) return { valid: false, operation: op, error: 'operation "molar_mass" requires temperatureK' }
        const solved = molarMassFromDensity(args.densityGPerL, args.pressureAtm, args.temperatureK)
        if (!solved.ok) return { valid: false, operation: op, error: solved.error }
        const out: GasLawResult = {
          valid: true,
          operation: op,
          pressureAtm: solved.pressureAtm,
          temperatureK: solved.temperatureK,
          temperatureC: solved.temperatureC,
          densityGPerL: solved.densityGPerL,
          molarMass: solved.molarMass,
          r: solved.r,
        }
        const mass = resolveMolarMass(args.formula, undefined)
        if (args.formula !== undefined && mass.formula !== null) out.formula = mass.formula
        return out
      }
      return { valid: false, error: `unknown operation "${op}" (expected ideal, molar_volume, density or molar_mass)` }
    },
  })

  const stoichiometry_tool = defineTool({
    name: 'stoichiometry',
    description: 'Limiting-reagent arithmetic for one balanced equation: give the equation and how much of each '
      + 'reactant you actually have (moles or grams), and get the extent of reaction, which reactant runs out '
      + 'first, how much of every product forms (mol and g) and what is left over. Optional `product` + '
      + '`actualYieldG` turn it into a theoretical / percent yield. Coefficients come from exact rational '
      + 'balancing, molar masses from IUPAC atomic weights — nothing is estimated. Every reactant in the '
      + 'equation must be given an amount: an unstated amount would make the limiting reagent undefined, so it '
      + 'is reported as an error instead of guessed.',
    parameters: {
      equation: { type: 'string', required: true, description: 'The equation (balanced or not), e.g. "C3H8 + O2 -> CO2 + H2O".' },
      reactants: {
        type: 'array',
        required: true,
        description: 'One entry per reactant in the equation, each with exactly one of moles or grams.',
        items: {
          type: 'object',
          additionalProperties: false,
          properties: {
            formula: { type: 'string', required: true, description: 'Reactant formula exactly as written in the equation, e.g. "C3H8".' },
            moles: { type: 'number', description: 'Amount available in mol (> 0).' },
            grams: { type: 'number', description: 'Mass available in g (> 0); converted with the formula molar mass.' },
          },
        },
      },
      product: { type: 'string', description: 'Product whose theoretical yield to report, e.g. "CO2".' },
      actualYieldG: { type: 'number', description: 'Measured mass of `product` in g — requires `product`; yields percentYield.' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          equation: { type: 'string' },
          extent: { type: 'number' },
          limiting: { type: 'array', items: { type: 'string' } },
          reactants: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                formula: { type: 'string', required: true },
                coefficient: { type: 'number', required: true },
                molarMass: { type: 'number', required: true },
                availableMoles: { type: 'number', required: true },
                availableGrams: { type: 'number', required: true },
                consumedMoles: { type: 'number', required: true },
                leftoverMoles: { type: 'number', required: true },
                leftoverGrams: { type: 'number', required: true },
                limiting: { type: 'boolean', required: true },
              },
            },
          },
          products: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                formula: { type: 'string', required: true },
                coefficient: { type: 'number', required: true },
                molarMass: { type: 'number', required: true },
                moles: { type: 'number', required: true },
                grams: { type: 'number', required: true },
              },
            },
          },
          yieldProduct: { type: 'string' },
          theoreticalYieldMol: { type: 'number' },
          theoreticalYieldG: { type: 'number' },
          actualYieldG: { type: 'number' },
          percentYield: { type: 'number' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderStoichiometry(value) }],
    },
    async execute(args: {
      equation: string
      reactants: { formula: string; moles?: number; grams?: number }[]
      product?: string
      actualYieldG?: number
    }): Promise<StoichiometryOut> {
      const result: StoichResult = stoichiometry({
        equation: args.equation,
        reactants: args.reactants ?? [],
        product: args.product,
        actualYieldG: args.actualYieldG,
      })
      if (!result.ok) return { valid: false, error: result.error }
      const out: StoichiometryOut = {
        valid: true,
        equation: result.equation,
        extent: result.extent,
        limiting: result.limiting,
        reactants: result.reactants,
        products: result.products,
      }
      if (result.yieldProduct !== undefined) out.yieldProduct = result.yieldProduct
      if (result.theoreticalYieldMol !== undefined) out.theoreticalYieldMol = result.theoreticalYieldMol
      if (result.theoreticalYieldG !== undefined) out.theoreticalYieldG = result.theoreticalYieldG
      if (result.actualYieldG !== undefined) out.actualYieldG = result.actualYieldG
      if (result.percentYield !== undefined) out.percentYield = result.percentYield
      return out
    },
  })

  const formula_from_composition_tool = defineTool({
    name: 'formula_from_composition',
    description: 'Empirical formula from an elemental analysis: give the mass percent (or the grams) of each '
      + 'element and get the empirical formula, the integer multiplier used, and — when `molarMass` is given — the '
      + 'molecular formula. This is the inverse of formula_mass and the arithmetic agents most often fumble: the '
      + 'smallest mole ratio is not necessarily 1 (Fe2O3 needs ×2 from 1 : 1.5), a measured molar mass is only a '
      + 'molecular formula if it is a whole multiple of the empirical mass, and a noisy analysis does not support '
      + 'integer subscripts at all. Subscripts come from the smallest multiplier k in [1, 8] whose scaled ratios '
      + 'are within ±0.05 of integers; when no k fits, the closest fit is returned with `exact: false` and a note '
      + 'instead of a quietly rounded answer. Formulas are printed in Hill order (C, H, then alphabetical).',
    parameters: {
      mode: {
        type: 'string',
        required: true,
        enum: ['percent', 'masses'],
        description: '"percent" = the values are mass percentages (they must sum to 100 ± 1); "masses" = grams.',
      },
      parts: {
        type: 'array',
        required: true,
        description: 'One entry per element: its symbol and the measured mass percent or mass in grams.',
        items: {
          type: 'object',
          additionalProperties: false,
          properties: {
            symbol: { type: 'string', required: true, description: 'Element symbol, e.g. "C", "Fe", "Na".' },
            value: { type: 'number', required: true, description: 'Mass percent (mode "percent") or grams (mode "masses"); must be > 0.' },
          },
        },
      },
      molarMass: { type: 'number', description: 'Measured molar mass in g/mol — when given, the molecular formula is reported too.' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          mode: { type: 'string' },
          empiricalFormula: { type: 'string' },
          empiricalMolarMass: { type: 'number' },
          multiplier: { type: 'number' },
          maxDeviation: { type: 'number' },
          exact: { type: 'boolean' },
          elements: {
            type: 'array',
            items: {
              type: 'object',
              additionalProperties: false,
              properties: {
                symbol: { type: 'string', required: true },
                massG: { type: 'number', required: true },
                atomicMass: { type: 'number', required: true },
                moles: { type: 'number', required: true },
                ratio: { type: 'number', required: true },
                subscript: { type: 'number', required: true },
                massFraction: { type: 'number', required: true },
              },
            },
          },
          molecularFormula: { type: 'string' },
          molecularMultiplier: { type: 'number' },
          molecularExact: { type: 'boolean' },
          measuredMolarMass: { type: 'number' },
          impliedMolarMass: { type: 'number' },
          molarMassErrorPercent: { type: 'number' },
          note: { type: 'string' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderComposition(value) }],
    },
    async execute(args: {
      mode: 'percent' | 'masses'
      parts: { symbol: string; value: number }[]
      molarMass?: number
    }): Promise<CompositionOut> {
      const result = compositionToFormula({ mode: args.mode, parts: args.parts ?? [], molarMass: args.molarMass })
      if (!result.ok) return { valid: false, error: result.error }
      const out: CompositionOut = {
        valid: true,
        mode: result.mode,
        empiricalFormula: result.empiricalFormula,
        empiricalMolarMass: result.empiricalMolarMass,
        multiplier: result.multiplier,
        maxDeviation: result.maxDeviation,
        exact: result.exact,
        elements: result.elements,
      }
      if (result.molecularFormula !== undefined) out.molecularFormula = result.molecularFormula
      if (result.molecularMultiplier !== undefined) out.molecularMultiplier = result.molecularMultiplier
      if (result.molecularExact !== undefined) out.molecularExact = result.molecularExact
      if (result.measuredMolarMass !== undefined) out.measuredMolarMass = result.measuredMolarMass
      if (result.impliedMolarMass !== undefined) out.impliedMolarMass = result.impliedMolarMass
      if (result.molarMassErrorPercent !== undefined) out.molarMassErrorPercent = result.molarMassErrorPercent
      if (result.note !== undefined) out.note = result.note
      return out
    },
  })

  const equilibrium_tool = defineTool({
    name: 'equilibrium',
    description: 'Weak acid / weak base equilibrium (ICE table) at 25 °C: give Ka (or Kb) and the formal '
      + 'concentration and get the amount dissociated (x, the positive root of x² + K·x − K·C = 0 — evaluated in '
      + 'the cancellation-free form 2KC/(K + √(K²+4KC))), [H+], [OH-], pH, pOH and the degree of dissociation. '
      + 'It also reports what the √(K·C) shortcut would have said and by how much it overstates x, plus a '
      + 'back-check K = x²/(C − x), so you can see when the approximation was not good enough (6.9 % off at '
      + '0.001 M ammonia with Kb = 1.8e-5). Monoprotic acids/bases only; no activity coefficients, ionic strength '
      + 'or polyprotic steps.',
    parameters: {
      mode: {
        type: 'string',
        required: true,
        enum: ['weak_acid', 'weak_base'],
        description: '"weak_acid" solves for [H+] (Ka); "weak_base" solves for [OH-] (Kb).',
      },
      k: { type: 'number', required: true, description: 'Ka (acid) or Kb (base), mol/L; must be > 0.' },
      concentration: {
        type: 'number',
        required: true,
        description: 'Formal (analytical) concentration of the acid/base in mol/L; must be > 0.',
      },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          mode: { type: 'string' },
          constantK: { type: 'number' },
          pK: { type: 'number' },
          concentration: { type: 'number' },
          x: { type: 'number' },
          hConc: { type: 'number' },
          ohConc: { type: 'number' },
          pH: { type: 'number' },
          pOH: { type: 'number' },
          degreeOfDissociation: { type: 'number' },
          percentDissociation: { type: 'number' },
          approximation: {
            type: 'object',
            additionalProperties: false,
            properties: {
              x: { type: 'number', required: true },
              errorPercent: { type: 'number', required: true },
            },
          },
          kFromResult: { type: 'number' },
          note: { type: 'string' },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderEquilibrium(value) }],
    },
    async execute(args: { mode: 'weak_acid' | 'weak_base'; k: number; concentration: number }): Promise<EquilibriumOut> {
      const result = weakEquilibrium({ mode: args.mode, k: args.k, concentration: args.concentration })
      if (!result.ok) return { valid: false, error: result.error }
      const out: EquilibriumOut = {
        valid: true,
        mode: result.mode,
        constantK: result.constantK,
        pK: result.pK,
        concentration: result.concentration,
        x: result.x,
        hConc: result.hConc,
        ohConc: result.ohConc,
        pH: result.pH,
        pOH: result.pOH,
        degreeOfDissociation: result.degreeOfDissociation,
        percentDissociation: result.percentDissociation,
        approximation: result.approximation,
        kFromResult: result.kFromResult,
      }
      if (result.note !== undefined) out.note = result.note
      return out
    },
  })

  const ksp_solubility_tool = defineTool({
    name: 'ksp_solubility',
    description: 'Solubility product (Ksp) → molar solubility, with the ion grouping made explicit. Give a salt '
      + 'formula and its Ksp and get the ions it breaks into (with counts and charges), the printed Ksp expression '
      + '(e.g. Ca3(PO4)2 → (3s)³(2s)² = 108·s⁵), molar solubility in mol/L and g/L, and the ion concentrations at '
      + 'saturation. Optionally pass a common ion already present (commonIon + commonIonConcentration) and it solves '
      + 'the shifted equilibrium numerically instead of pretending the answer is unchanged. The grouping is verified, '
      + 'not guessed: a formula that matches two different ion splits (Hg2Cl2 is both Hg2²⁺·2Cl⁻ and 2Hg⁺·2Cl⁻) is '
      + 'reported as ambiguous with the candidates, and hydrates or charged formulas are rejected with the reason.',
    parameters: {
      formula: { type: 'string', required: true, description: 'Neutral salt formula, e.g. "AgCl", "Ca(OH)2", "Ca3(PO4)2".' },
      ksp: { type: 'number', required: true, description: 'Solubility product at the temperature it was measured, mol/L units as tabulated; must be > 0.' },
      cation: { type: 'string', description: 'Force the cation when the split is ambiguous, e.g. "Hg^2+", "Ca2+", "NH4+".' },
      anion: { type: 'string', description: 'Force the anion when the split is ambiguous, e.g. "Cl-", "PO4^3-".' },
      commonIon: { type: 'string', description: 'One of the salt\'s own ions already present, e.g. "Cl-" for AgCl in brine.' },
      commonIonConcentration: { type: 'number', description: 'Concentration of that common ion in mol/L; required with commonIon.' },
    },
    output: {
      schema: {
        type: 'object',
        additionalProperties: false,
        properties: {
          valid: { type: 'boolean', required: true },
          formula: { type: 'string' },
          ksp: { type: 'number' },
          pKsp: { type: 'number' },
          cation: {
            type: 'object',
            additionalProperties: false,
            properties: {
              symbol: { type: 'string', required: true },
              charge: { type: 'number', required: true },
              count: { type: 'number', required: true },
              concentration: { type: 'number', required: true },
            },
          },
          anion: {
            type: 'object',
            additionalProperties: false,
            properties: {
              symbol: { type: 'string', required: true },
              charge: { type: 'number', required: true },
              count: { type: 'number', required: true },
              concentration: { type: 'number', required: true },
            },
          },
          stoichiometricFactor: { type: 'number' },
          ionCount: { type: 'number' },
          expression: { type: 'string' },
          molarMass: { type: 'number' },
          molarSolubility: { type: 'number' },
          massSolubility: { type: 'number' },
          chargeBalance: { type: 'number' },
          commonIon: {
            type: 'object',
            additionalProperties: false,
            properties: {
              symbol: { type: 'string', required: true },
              concentration: { type: 'number', required: true },
              molarSolubility: { type: 'number', required: true },
              massSolubility: { type: 'number', required: true },
              reductionFactor: { type: 'number', required: true },
              residual: { type: 'number', required: true },
            },
          },
          notes: { type: 'array', items: { type: 'string' } },
          candidates: { type: 'array', items: { type: 'string' } },
          error: { type: 'string' },
        },
      },
      render: (_args: unknown, value: unknown) => [{ type: 'text', text: renderKsp(value) }],
    },
    async execute(args: {
      formula: string
      ksp: number
      cation?: string
      anion?: string
      commonIon?: string
      commonIonConcentration?: number
    }): Promise<KspOut> {
      const result = kspSolubility({
        formula: args.formula,
        ksp: args.ksp,
        cation: args.cation,
        anion: args.anion,
        commonIon: args.commonIon,
        commonIonConcentration: args.commonIonConcentration,
      })
      if (!result.ok) {
        const out: KspOut = { valid: false, error: result.error }
        if (result.candidates !== undefined) out.candidates = result.candidates
        return out
      }
      const out: KspOut = {
        valid: true,
        formula: result.formula,
        ksp: result.ksp,
        pKsp: result.pKsp,
        cation: result.cation,
        anion: result.anion,
        stoichiometricFactor: result.stoichiometricFactor,
        ionCount: result.ionCount,
        expression: result.expression,
        molarMass: result.molarMass,
        molarSolubility: result.molarSolubility,
        massSolubility: result.massSolubility,
        chargeBalance: result.chargeBalance,
        notes: result.notes,
      }
      if (result.commonIon !== undefined) out.commonIon = result.commonIon
      return out
    },
  })

  return {
    formula_mass,
    mass_moles,
    solution,
    ph,
    balance_equation,
    gas_law,
    stoichiometry: stoichiometry_tool,
    formula_from_composition: formula_from_composition_tool,
    equilibrium: equilibrium_tool,
    ksp_solubility: ksp_solubility_tool,
  }
}
