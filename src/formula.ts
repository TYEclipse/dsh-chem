/**
 * Composition → formula: turn the elemental analysis of a compound (mass
 * percentages, or plain masses) into an empirical formula, and — when a molar
 * mass is known — the molecular formula.
 *
 * The arithmetic is the classic divide-by-atomic-mass / divide-by-smallest /
 * scale-to-integers loop, which agents get wrong in predictable ways: they
 * forget that the smallest ratio itself may not be 1, they round 1.5 to 2
 * without noticing, and they happily report a "molecular formula" from a molar
 * mass that is not an integer multiple of the empirical mass. This module
 * reports the multiplier, the residual and an explicit `exact` flag instead of
 * quietly rounding.
 *
 * Deterministic rules (documented in the README):
 *   - `percent` is mass percent and must sum to 100 ± 1 (a typo is an error,
 *     not something to normalise away silently);
 *   - subscripts come from the SMALLEST multiplier k ∈ [1, 8] whose scaled
 *     ratios are all within ±0.05 of an integer;
 *   - if no k fits, the closest fit is reported with `exact: false`;
 *   - formulas are written in Hill order (C, H, then alphabetical; alphabetical
 *     when there is no carbon) so the same analysis always prints the same
 *     string regardless of input order.
 *
 * @module dsh-chem/formula
 */

import { formulaMass, parseFormula, roundOut } from './core.ts'
import { elementBySymbolLenient } from './periodic.ts'

/** Largest number of elements accepted in one composition. */
export const MAX_ELEMENTS = 12

/** Largest integer multiplier searched for (empirical subscripts are coprime). */
export const MAX_MULTIPLIER = 8

/** |ratio × k − round(ratio × k)| must be within this to count as an integer. */
export const RATIO_TOLERANCE = 0.05

/** A molar mass within this relative distance of an integer multiple is accepted. */
export const MOLECULAR_TOLERANCE = 0.02

export interface CompositionPart {
  /** Element symbol, e.g. "C" (single-letter lowercase is corrected, as elsewhere). */
  symbol: string
  /** Mass percent (mode "percent") or mass in grams (mode "masses"). */
  value: number
}

export interface CompositionInput {
  /** "percent" = the values are mass percentages (must sum to 100 ± 1); "masses" = grams. */
  mode: 'percent' | 'masses'
  parts: CompositionPart[]
  /** Measured molar mass in g/mol — turns the answer into a molecular formula. */
  molarMass?: number
}

export interface CompositionElementOut {
  symbol: string
  /** Mass used in the calculation (grams, or grams per 100 g for `percent`). */
  massG: number
  atomicMass: number
  moles: number
  /** moles ÷ smallest moles (before scaling). */
  ratio: number
  /** Integer subscript in the empirical formula. */
  subscript: number
  /** Mass percent of this element in the empirical formula (round-trip check). */
  massFraction: number
}

export type CompositionResult =
  | {
      ok: true
      mode: 'percent' | 'masses'
      empiricalFormula: string
      empiricalMolarMass: number
      multiplier: number
      maxDeviation: number
      exact: boolean
      elements: CompositionElementOut[]
      molecularFormula?: string
      molecularMultiplier?: number
      molecularExact?: boolean
      measuredMolarMass?: number
      impliedMolarMass?: number
      molarMassErrorPercent?: number
      note?: string
    }
  | { ok: false; error: string }

interface ElementMass {
  symbol: string
  atomicMass: number
  mass: number
}

/** Hill-system ordering: C, then H, then alphabetical; alphabetical without C. */
export function hillOrder(symbols: string[]): string[] {
  const unique = [...new Set(symbols)]
  if (unique.includes('C')) {
    const rest = unique.filter((s) => s !== 'C' && s !== 'H').sort()
    return ['C', ...(unique.includes('H') ? ['H'] : []), ...rest]
  }
  return unique.sort()
}

/** Render subscripts as a formula string in Hill order (subscript 1 omitted). */
export function renderFormula(counts: Record<string, number>): string {
  return hillOrder(Object.keys(counts))
    .map((symbol) => {
      const count = counts[symbol]!
      return count === 1 ? symbol : `${symbol}${count}`
    })
    .join('')
}

/** Molar mass of a rendered formula, via the plugin's own IUPAC table. */
function massOfFormula(formula: string): number | { error: string } {
  const parsed = parseFormula(formula)
  if (!parsed.ok) return { error: parsed.error }
  return formulaMass(parsed).molarMass
}

/**
 * Resolve an elemental analysis → empirical (and molecular) formula.
 *
 * `mode` says what the numbers mean ("percent" = mass %, "masses" = grams);
 * `molarMass` is optional.
 */
export function compositionToFormula(input: CompositionInput): CompositionResult {
  const mode = input.mode
  if (mode !== 'percent' && mode !== 'masses') {
    return { ok: false, error: `unknown mode "${mode}" (expected percent or masses)` }
  }
  const parts = input.parts
  if (!Array.isArray(parts)) return { ok: false, error: '`parts` must be an array of { symbol, value }' }
  if (parts.length === 0) return { ok: false, error: '`parts` must name at least one element' }
  if (parts.length > MAX_ELEMENTS) {
    return { ok: false, error: `too many elements (${parts.length}); at most ${MAX_ELEMENTS} per analysis` }
  }

  const resolved = new Map<string, ElementMass>()
  let total = 0
  for (const part of parts) {
    const raw = part?.symbol
    const value = part?.value
    if (typeof raw !== 'string' || raw.trim() === '') {
      return { ok: false, error: `every part needs a symbol (got ${JSON.stringify(raw)})` }
    }
    if (typeof value !== 'number' || !Number.isFinite(value)) {
      return { ok: false, error: `"${raw}": value must be a finite number (got ${value})` }
    }
    if (value <= 0) return { ok: false, error: `"${raw}": value must be positive (got ${value})` }
    const element = elementBySymbolLenient(raw)
    if (!element) return { ok: false, error: `unknown element symbol "${raw}"` }
    if (resolved.has(element.symbol)) {
      return { ok: false, error: `element ${element.symbol} given twice ("${raw}"); merge the entries first` }
    }
    resolved.set(element.symbol, { symbol: element.symbol, atomicMass: element.atomicMass, mass: value })
    total += value
  }

  if (mode === 'percent' && Math.abs(total - 100) > 1) {
    return {
      ok: false,
      error: `mass percentages sum to ${roundOut(total, 6)}, expected 100 ± 1 — check for a missing or duplicated element`,
    }
  }

  const rows = [...resolved.values()]
  const moles = rows.map((row) => row.mass / row.atomicMass)
  const smallest = Math.min(...moles)
  const ratios = moles.map((value) => value / smallest)

  let multiplier = 1
  let bestDeviation = Number.POSITIVE_INFINITY
  let exact = false
  for (let k = 1; k <= MAX_MULTIPLIER; k += 1) {
    const deviation = Math.max(...ratios.map((ratio) => Math.abs(ratio * k - Math.round(ratio * k))))
    if (deviation < bestDeviation - 1e-15) {
      multiplier = k
      bestDeviation = deviation
    }
    if (deviation <= RATIO_TOLERANCE) {
      multiplier = k
      bestDeviation = deviation
      exact = true
      break
    }
  }

  const counts: Record<string, number> = {}
  const subscripts = ratios.map((ratio) => Math.max(1, Math.round(ratio * multiplier)))
  rows.forEach((row, index) => {
    counts[row.symbol] = subscripts[index]!
  })
  const empiricalFormula = renderFormula(counts)
  const empiricalMass = massOfFormula(empiricalFormula)
  if (typeof empiricalMass === 'object') {
    return { ok: false, error: `could not re-read the empirical formula "${empiricalFormula}": ${empiricalMass.error}` }
  }

  const elements: CompositionElementOut[] = rows.map((row, index) => {
    const subscript = subscripts[index]!
    const contribution = subscript * row.atomicMass
    return {
      symbol: row.symbol,
      massG: roundOut(row.mass),
      atomicMass: roundOut(row.atomicMass),
      moles: roundOut(moles[index]!),
      ratio: roundOut(ratios[index]!),
      subscript,
      massFraction: roundOut((contribution / empiricalMass) * 100),
    }
  })
  const order = hillOrder(Object.keys(counts))
  elements.sort((a, b) => order.indexOf(a.symbol) - order.indexOf(b.symbol))

  const result: CompositionResult = {
    ok: true,
    mode,
    empiricalFormula,
    empiricalMolarMass: roundOut(empiricalMass),
    multiplier,
    maxDeviation: roundOut(bestDeviation),
    exact,
    elements,
  }

  if (!exact) {
    result.note = `no integer ratio found within ±${RATIO_TOLERANCE} (closest fit at k = ${multiplier}, `
      + `max deviation ${roundOut(bestDeviation, 6)}) — the subscripts shown are the closest fit, not a formula `
      + 'that the data actually supports; re-check the analysis or the element list'
  }

  if (input.molarMass !== undefined) {
    if (!Number.isFinite(input.molarMass) || input.molarMass <= 0) {
      return { ok: false, error: `molarMass must be positive (got ${input.molarMass})` }
    }
    const multiple = input.molarMass / empiricalMass
    const nearest = Math.max(1, Math.round(multiple))
    const relative = Math.abs(multiple - nearest) / nearest
    result.measuredMolarMass = roundOut(input.molarMass)
    result.molecularMultiplier = nearest
    result.molecularExact = relative <= MOLECULAR_TOLERANCE
    result.impliedMolarMass = roundOut(nearest * empiricalMass)
    result.molarMassErrorPercent = roundOut((Math.abs(nearest * empiricalMass - input.molarMass) / input.molarMass) * 100)
    const molecularCounts: Record<string, number> = {}
    for (const [symbol, count] of Object.entries(counts)) molecularCounts[symbol] = count * nearest
    result.molecularFormula = renderFormula(molecularCounts)
    if (!result.molecularExact) {
      const extra = `${input.molarMass} g/mol is ${roundOut(multiple, 6)} × the empirical mass `
        + `(${roundOut(empiricalMass, 6)} g/mol), which is not close to a whole number within `
        + `${MOLECULAR_TOLERANCE * 100} % — "${result.molecularFormula}" is the nearest whole multiple, not a result `
        + 'of your data'
      result.note = result.note ? `${result.note}; ${extra}` : extra
    }
  }

  return result
}
