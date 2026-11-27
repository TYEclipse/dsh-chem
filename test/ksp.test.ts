/**
 * Tests for `ksp_solubility` (solubility product → molar solubility).
 *
 * Every numeric anchor comes from the independent Python oracle
 * `test/anchors-chem.py` (its own ion table, matching salts by element counts
 * instead of splitting the parsed groups, and Newton instead of bisection for
 * the common-ion case) — see the `ksp.*` lines in its output.
 *
 * ORACLE: test/anchors-chem.py
 *
 * @module dsh-chem/test/ksp
 */

import { describe, expect, it } from 'vitest'
import { buildChemTools } from '../src/tools.ts'
import { assertLossless } from './helpers/lossless-json.ts'

const tools = buildChemTools()
const executeKsp = tools.ksp_solubility.execute as (args: Record<string, unknown>) => Promise<Record<string, unknown>>

/** Build args with only defined keys (args-side lossless JSON discipline). */
function pick(entries: Record<string, unknown>): Record<string, unknown> {
  const result: Record<string, unknown> = {}
  for (const [key, value] of Object.entries(entries)) {
    if (value !== undefined) result[key] = value
  }
  return result
}

function relErr(actual: number, expected: number): number {
  return Math.abs(actual - expected) / Math.abs(expected)
}

interface IonOut { symbol: string; charge: number; count: number; concentration: number }

/* ------------------------------------------------------------------ */
/* groupings + molar solubility (anchors from the python oracle)       */
/* ------------------------------------------------------------------ */

const SOLUBILITY_ANCHORS: Array<{
  salt: string
  ksp: number
  cation: string
  cationCount: number
  anion: string
  anionCount: number
  factor: number
  order: number
  s: number
  massSolubility: number
  molarMass: number
}> = [
  { salt: 'AgCl', ksp: 1.8e-10, cation: 'Ag^+', cationCount: 1, anion: 'Cl^-', anionCount: 1, factor: 1, order: 2, s: 1.3416407864998738e-05, massSolubility: 0.001922839575211619, molarMass: 143.32 },
  { salt: 'Ca(OH)2', ksp: 5.5e-6, cation: 'Ca^2+', cationCount: 1, anion: 'OH^-', anionCount: 2, factor: 4, order: 3, s: 0.01111990045284658, massSolubility: 0.8238956643523088, molarMass: 74.092 },
  { salt: 'PbCl2', ksp: 1.7e-5, cation: 'Pb^2+', cationCount: 1, anion: 'Cl^-', anionCount: 2, factor: 4, order: 3, s: 0.01619805900638742, massSolubility: 4.504680209676342, molarMass: 278.1 },
  { salt: 'Ag2CrO4', ksp: 1.1e-12, cation: 'Ag^+', cationCount: 2, anion: 'CrO4^2-', anionCount: 1, factor: 4, order: 3, s: 6.502957234256938e-05, massSolubility: 0.021572390092345226, molarMass: 331.732 },
  { salt: 'BaSO4', ksp: 1.1e-10, cation: 'Ba^2+', cationCount: 1, anion: 'SO4^2-', anionCount: 1, factor: 1, order: 2, s: 1.0488088481701516e-05, massSolubility: 0.00244777301839039, molarMass: 233.386 },
  { salt: 'Mg(OH)2', ksp: 5.6e-12, cation: 'Mg^2+', cationCount: 1, anion: 'OH^-', anionCount: 2, factor: 4, order: 3, s: 0.00011186889420813974, massSolubility: 0.006524082041324502, molarMass: 58.319 },
  { salt: 'Ca3(PO4)2', ksp: 2.07e-33, cation: 'Ca^2+', cationCount: 3, anion: 'PO4^3-', anionCount: 2, factor: 108, order: 5, s: 1.138962218587898e-07, massSolubility: 3.532764671882827e-05, molarMass: 310.174 },
  { salt: 'Fe(OH)3', ksp: 2.79e-39, cation: 'Fe^3+', cationCount: 1, anion: 'OH^-', anionCount: 3, factor: 27, order: 4, s: 1.0082311468433848e-10, massSolubility: 1.0774562973856517e-08, molarMass: 106.866 },
  { salt: 'LiF', ksp: 1.8e-3, cation: 'Li^+', cationCount: 1, anion: 'F^-', anionCount: 1, factor: 1, order: 2, s: 0.042426406871192854, massSolubility: 1.1004561414250003, molarMass: 25.938 },
  { salt: 'NH4NO3', ksp: 1.0e-3, cation: 'NH4^+', cationCount: 1, anion: 'NO3^-', anionCount: 1, factor: 1, order: 2, s: 0.03162277660168379, massSolubility: 2.531181907528576, molarMass: 80.043 },
]

describe('ksp_solubility: ion grouping and molar solubility', () => {
  it('splits each salt into the ions the oracle found and matches its solubility', async () => {
    for (const anchor of SOLUBILITY_ANCHORS) {
      const result = await executeKsp(pick({ formula: anchor.salt, ksp: anchor.ksp }))
      assertLossless(result)
      expect(result.valid, anchor.salt).toBe(true)
      const cation = result.cation as IonOut
      const anion = result.anion as IonOut
      expect(cation.symbol, anchor.salt).toBe(anchor.cation)
      expect(cation.count, anchor.salt).toBe(anchor.cationCount)
      expect(anion.symbol, anchor.salt).toBe(anchor.anion)
      expect(anion.count, anchor.salt).toBe(anchor.anionCount)
      expect(result.stoichiometricFactor, anchor.salt).toBe(anchor.factor)
      expect(result.ionCount, anchor.salt).toBe(anchor.order)
      expect(relErr(result.molarSolubility as number, anchor.s)).toBeLessThan(1e-12)
      expect(relErr(result.massSolubility as number, anchor.massSolubility)).toBeLessThan(1e-9)
      expect(relErr(result.molarMass as number, anchor.molarMass)).toBeLessThan(1e-9)
      // the concentrations are the stoichiometric multiples of s
      expect(relErr(cation.concentration, anchor.cationCount * anchor.s)).toBeLessThan(1e-12)
      expect(relErr(anion.concentration, anchor.anionCount * anchor.s)).toBeLessThan(1e-12)
      // charge always balances — that is what makes the grouping a salt at all
      expect(result.chargeBalance, anchor.salt).toBe(0)
    }
  })

  it('prints the Ksp expression in the form a reader can check by hand', async () => {
    const cases: Array<{ salt: string; ksp: number; expression: string }> = [
      { salt: 'AgCl', ksp: 1.8e-10, expression: 'Ksp = [Ag^+]·[Cl^-] = s·s = s^2' },
      { salt: 'Ca(OH)2', ksp: 5.5e-6, expression: 'Ksp = [Ca^2+]·[OH^-]^2 = s·(2s) = 4·s^3' },
      { salt: 'Ca3(PO4)2', ksp: 2.07e-33, expression: 'Ksp = [Ca^2+]^3·[PO4^3-]^2 = (3s)·(2s) = 108·s^5' },
      { salt: 'BaSO4', ksp: 1.1e-10, expression: 'Ksp = [Ba^2+]·[SO4^2-] = s·s = s^2' },
    ]
    for (const item of cases) {
      const result = await executeKsp(pick({ formula: item.salt, ksp: item.ksp }))
      expect(result.expression, item.salt).toBe(item.expression)
    }
  })

  it('reports pKsp and a note when the answer is too concentrated for ideal-solution Ksp', async () => {
    const dilute = await executeKsp(pick({ formula: 'AgCl', ksp: 1.8e-10 }))
    expect(relErr(dilute.pKsp as number, 9.744727494896693)).toBeLessThan(1e-9)
    const notes = dilute.notes as string[]
    expect(notes.some((n) => n.includes('temperature dependent'))).toBe(true)
    expect(notes.some((n) => n.includes('activity coefficients'))).toBe(false)

    // s = 1 mol/L here: the ideal-solution assumption is the thing to flag
    const concentrated = await executeKsp(pick({ formula: 'AgCl', ksp: 1 }))
    const loudNotes = concentrated.notes as string[]
    expect(loudNotes.some((n) => n.includes('activity coefficients'))).toBe(true)
  })

  it('does not report a common-ion block when none was asked for', async () => {
    const result = await executeKsp(pick({ formula: 'AgCl', ksp: 1.8e-10 }))
    expect('commonIon' in result).toBe(false)
    assertLossless(result)
  })
})

/* ------------------------------------------------------------------ */
/* ambiguity + honest refusals                                          */
/* ------------------------------------------------------------------ */

describe('ksp_solubility: ambiguous groupings are reported, not guessed', () => {
  it('refuses Hg2Cl2 until the mercury(I) form is chosen', async () => {
    const ambiguous = await executeKsp(pick({ formula: 'Hg2Cl2', ksp: 1.3e-18 }))
    assertLossless(ambiguous)
    expect(ambiguous.valid).toBe(false)
    expect(String(ambiguous.error)).toContain('matches 2 ion groupings')
    const candidates = ambiguous.candidates as string[]
    expect(candidates).toContain('1 Hg2^2+ + 2 Cl^-')
    expect(candidates).toContain('2 Hg^+ + 2 Cl^-')
    // no answer fields leak out of a refusal
    expect(candidates.length).toBe(2)
    expect('molarSolubility' in ambiguous).toBe(false)
  })

  it('accepts the grouped mercury(I) ion once given', async () => {
    const result = await executeKsp(pick({ formula: 'Hg2Cl2', ksp: 1.3e-18, cation: 'Hg2^2+' }))
    assertLossless(result)
    expect(result.valid).toBe(true)
    expect((result.cation as IonOut).symbol).toBe('Hg2^2+')
    expect((result.anion as IonOut).count).toBe(2)
    expect(relErr(result.molarSolubility as number, 6.87534433537071e-07)).toBeLessThan(1e-12)
    expect(relErr(result.massSolubility as number, 0.0003245712553841805)).toBeLessThan(1e-9)
  })

  it('still refuses when only the anion is pinned (both readings use Cl^-)', async () => {
    const result = await executeKsp(pick({ formula: 'Hg2Cl2', ksp: 1.3e-18, anion: 'Cl^-' }))
    expect(result.valid).toBe(false)
    expect(String(result.error)).toContain('give both cation and anion')
  })
})

describe('ksp_solubility: rejections name the reason', () => {
  const rejects: Array<{ label: string; args: Record<string, unknown>; needle: string }> = [
    { label: 'zero Ksp', args: { formula: 'AgCl', ksp: 0 }, needle: 'ksp must be > 0' },
    { label: 'negative Ksp', args: { formula: 'AgCl', ksp: -1e-10 }, needle: 'ksp must be > 0' },
    { label: 'empty formula', args: { formula: '  ', ksp: 1e-10 }, needle: 'formula is empty' },
    { label: 'hydrate', args: { formula: 'CuSO4·5H2O', ksp: 1e-10 }, needle: 'is a hydrate' },
    { label: 'leading coefficient', args: { formula: '2NaCl', ksp: 1e-10 }, needle: 'starts with a coefficient' },
    { label: 'charged formula', args: { formula: 'Ag+', ksp: 1e-10 }, needle: 'carries a charge' },
    { label: 'single group', args: { formula: 'Na', ksp: 1e-10 }, needle: 'a salt needs a cation and an anion' },
    { label: 'mixed-valence oxide', args: { formula: 'Fe3O4', ksp: 1e-10 }, needle: 'does not split into two recognised ions' },
    { label: 'unparseable', args: { formula: 'Xx2Cl', ksp: 1e-10 }, needle: 'unknown element symbol' },
    { label: 'unknown forced cation', args: { formula: 'AgCl', ksp: 1e-10, cation: 'Zz^+' }, needle: 'is not in the recognised ion table' },
    { label: 'forced cation from another salt', args: { formula: 'AgCl', ksp: 1e-10, cation: 'K^+' }, needle: 'no ion grouping' },
    { label: 'common ion of another salt', args: { formula: 'AgCl', ksp: 1e-10, commonIon: 'SO4^2-', commonIonConcentration: 0.1 }, needle: 'is not one of the ions' },
    { label: 'unknown common ion', args: { formula: 'AgCl', ksp: 1e-10, commonIon: 'Zz^+', commonIonConcentration: 0.1 }, needle: 'is not in the recognised ion table' },
    { label: 'common ion without a concentration', args: { formula: 'AgCl', ksp: 1e-10, commonIon: 'Cl^-' }, needle: 'commonIonConcentration must be a finite number' },
    { label: 'common ion at zero concentration', args: { formula: 'AgCl', ksp: 1e-10, commonIon: 'Cl^-', commonIonConcentration: 0 }, needle: 'must be > 0' },
  ]

  it.each(rejects)('$label', async ({ args, needle }) => {
    const result = await executeKsp(pick(args))
    assertLossless(result)
    expect(result.valid).toBe(false)
    expect(String(result.error)).toContain(needle)
  })

  it('never leaks numeric fields on a refusal', async () => {
    const result = await executeKsp(pick({ formula: 'Fe3O4', ksp: 1e-10 }))
    for (const key of ['molarSolubility', 'massSolubility', 'stoichiometricFactor', 'cation', 'anion', 'pKsp']) {
      expect(key in result, key).toBe(false)
    }
  })

  it('rejects a non-finite Ksp at the schema layer (before the maths runs)', async () => {
    await expect(executeKsp(pick({ formula: 'AgCl', ksp: Number.NaN }))).rejects.toThrow(/finite/)
  })
})

/* ------------------------------------------------------------------ */
/* common-ion effect (anchors from the python oracle's Newton solve)    */
/* ------------------------------------------------------------------ */

describe('ksp_solubility: common-ion effect', () => {
  const anchors: Array<{ common: string; c: number; s: number; reduction: number; massSolubility: number }> = [
    { common: 'Cl^-', c: 0.01, s: 1.7999967600116638e-08, reduction: 745.3573341383016, massSolubility: 2.5797553564487162e-06 },
    { common: 'Ag^+', c: 0.01, s: 1.7999967600116638e-08, reduction: 745.3573341383016, massSolubility: 2.5797553564487162e-06 },
  ]

  it.each(anchors)('AgCl in 0.01 M $common', async (anchor) => {
    const result = await executeKsp(pick({
      formula: 'AgCl',
      ksp: 1.8e-10,
      commonIon: anchor.common,
      commonIonConcentration: anchor.c,
    }))
    assertLossless(result)
    expect(result.valid).toBe(true)
    const dense = result.commonIon as { symbol: string; concentration: number; molarSolubility: number; massSolubility: number; reductionFactor: number; residual: number }
    expect(dense.symbol).toBe(anchor.common)
    expect(dense.concentration).toBe(anchor.c)
    expect(relErr(dense.molarSolubility, anchor.s)).toBeLessThan(1e-9)
    expect(relErr(dense.massSolubility, anchor.massSolubility)).toBeLessThan(1e-9)
    expect(relErr(dense.reductionFactor, anchor.reduction)).toBeLessThan(1e-6)
    // the solver's own back-substitution residual must be reported and tiny
    expect(Math.abs(dense.residual)).toBeLessThan(1e-9)
    // and the plain-water answer is still reported next to it
    expect(relErr(result.molarSolubility as number, 1.3416407864998738e-05)).toBeLessThan(1e-12)
  })

  it('solves the shifted equilibrium exactly, not by dividing Ksp by the common ion', async () => {
    const c = 0.01
    const result = await executeKsp(pick({ formula: 'AgCl', ksp: 1.8e-10, commonIon: 'Cl^-', commonIonConcentration: c }))
    const dense = result.commonIon as { molarSolubility: number }
    const s = dense.molarSolubility
    // (s)(s + C) = Ksp must hold to the last few bits
    expect(relErr(s * (s + c), 1.8e-10)).toBeLessThan(1e-12)
    // the naive Ksp/C is a different (2e-4 relative) answer — the tool must not use it
    expect(Math.abs(s - 1.8e-8) / 1.8e-8).toBeGreaterThan(1e-7)
  })

  it('handles a common ion for a non-1:1 salt', async () => {
    // PbCl2 with 0.1 M Cl-: (s)(2s + 0.1)^2 = 1.7e-5
    const result = await executeKsp(pick({ formula: 'PbCl2', ksp: 1.7e-5, commonIon: 'Cl^-', commonIonConcentration: 0.1 }))
    const dense = result.commonIon as { molarSolubility: number; residual: number; reductionFactor: number }
    const s = dense.molarSolubility
    expect(relErr(s * (2 * s + 0.1) ** 2, 1.7e-5)).toBeLessThan(1e-9)
    expect(dense.reductionFactor).toBeGreaterThan(1)
    expect(Math.abs(dense.residual)).toBeLessThan(1e-9)
  })
})
