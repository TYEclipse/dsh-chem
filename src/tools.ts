/**
 * Tool definitions for dsh-chem: four deterministic chemistry tools exposed to
 * every agent via defineTool. Strict JSON-schema parameter surfaces, explicit
 * result interfaces matching the inferred output types.
 *
 * @module dsh-chem/tools
 */

import { defineTool, type ToolDefinition } from '@deepseek-ai/dsh-tools'
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

  return { formula_mass, mass_moles, solution, ph }
}
