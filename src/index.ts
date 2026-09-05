/**
 * dsh-chem — deterministic chemistry math for DeepSeek Harness.
 *
 * Four zero-dependency tools (pure arithmetic over IUPAC data):
 *   formula_mass — parse a chemical formula and compute its molar mass + breakdown
 *   mass_moles   — convert mass (g) <-> amount (mol) <-> particle count
 *   solution     — molarity, dilution (C1V1 = C2V2) and percent solutions
 *   ph           — pH / pOH / [H+] / [OH-] conversions + Henderson-Hasselbalch buffers
 *
 * Agents get chemistry arithmetic wrong in predictable ways: molar-mass sums,
 * case-sensitivity of element symbols (Co vs CO), dilution factors and log/antilog
 * pH math. These tools replace all of it with deterministic tables and formulas.
 *
 * @module dsh-chem
 */

import type { Context } from '@deepseek-ai/cordis'
import type { Agent } from '@deepseek-ai/dsh-agent'
import z from '@deepseek-ai/schemastery'
import { buildChemTools, type ToolSet } from './tools.ts'

/** Stable Cordis plugin name (also the config key under `plugins:`). */
export const name = 'dsh-chem'

/** Services required before tool registration can start. */
export const inject = ['agents', 'tools']

/** Plugin configuration (no tunables today — kept for schema completeness). */
export interface Config {
  /** Reserved for future options. */
  _?: never
}

export const Config: z<Config> = z.object({})

/** Register every chemistry tool on one agent; returns the disposer. */
function decorate(agent: Agent, tools: ToolSet): () => void {
  const disposers = Object.values(tools).map((definition) => agent.ctx.tools.register(definition))
  return () => {
    for (const dispose of disposers) {
      try {
        dispose()
      } catch {
        // already disposed
      }
    }
  }
}

/** Mount the chemistry tools on every live agent and every future one. */
export function apply(ctx: Context, _config: Config): void {
  const tools = buildChemTools()
  const disposers = new Set<() => void>()

  const decorateAgent = (agent: Agent): void => {
    try {
      disposers.add(decorate(agent, tools))
    } catch (error) {
      ctx.logger('chem').warn(`tool registration for agent ${agent.id} failed: ${error instanceof Error ? error.message : String(error)}`)
    }
  }

  for (const agent of ctx.agents.list()) decorateAgent(agent)
  const off = ctx.on('agent/created', ({ agent }) => decorateAgent(agent))

  ctx.effect(() => () => {
    off()
    for (const dispose of disposers) {
      try {
        dispose()
      } catch {
        // already disposed
      }
    }
    disposers.clear()
  })
}
