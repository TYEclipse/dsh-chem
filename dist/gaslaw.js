/**
 * Gas-law core for dsh-chem: ideal-gas arithmetic with unit conversions.
 *
 * Fixed unit surface (documented, no unit guessing):
 *   pressure    atm   (also reported in kPa and mmHg)
 *   volume      L     (also reported in mL)
 *   amount      mol
 *   temperature K     (also reported in °C)
 *
 * Molar gas constant: R = 8.31446261815324 J/(mol·K) (CODATA 2018) and
 * 1 atm = 101.325 kPa, therefore R = 8.31446261815324 / 101.325 L·atm/(mol·K).
 * Cross-validated against the published molar volume of an ideal gas at STP
 * (0 °C, 1 atm) = 22.414 L/mol.
 *
 * @module dsh-chem/gaslaw
 */
/** Molar gas constant in J/(mol·K) — CODATA 2018 exact value. */
export const R_JOULES = 8.31446261815324;
/** 1 atm in kPa (exact by definition of the standard atmosphere). */
export const ATM_TO_KPA = 101.325;
/** 1 atm in mmHg (conventional). */
export const ATM_TO_MMHG = 760;
/** Molar gas constant in L·atm/(mol·K). */
export const R_ATM_L = R_JOULES / ATM_TO_KPA;
/** Kelvin offset from degrees Celsius. */
export const KELVIN_OFFSET = 273.15;
function positive(value, label) {
    if (value === undefined || !Number.isFinite(value))
        return { ok: false, error: `${label} is required and must be a finite number` };
    if (value <= 0)
        return { ok: false, error: `${label} must be positive (got ${value})` };
    return value;
}
function withConversions(pressureAtm, volumeL, moles, temperatureK) {
    return {
        pressureAtm,
        pressureKpa: pressureAtm * ATM_TO_KPA,
        pressureMmhg: pressureAtm * ATM_TO_MMHG,
        volumeL,
        volumeMl: volumeL * 1000,
        moles,
        temperatureK,
        temperatureC: temperatureK - KELVIN_OFFSET,
        r: R_ATM_L,
    };
}
/** Solve PV = nRT for whichever of the four quantities was left out. */
export function idealGas(input) {
    const given = [];
    for (const key of ['pressureAtm', 'volumeL', 'moles', 'temperatureK']) {
        if (input[key] !== undefined) {
            given.push(key === 'pressureAtm' ? 'pressure' : key === 'volumeL' ? 'volume' : key === 'moles' ? 'moles' : 'temperature');
        }
    }
    if (given.length !== 3) {
        return {
            ok: false,
            error: `exactly three of pressureAtm, volumeL, moles and temperatureK are required (got ${given.length}: ${given.join(', ') || 'none'})`,
        };
    }
    const checks = [
        ['pressureAtm', input.pressureAtm],
        ['volumeL', input.volumeL],
        ['moles', input.moles],
        ['temperatureK', input.temperatureK],
    ];
    for (const [label, value] of checks) {
        if (value === undefined)
            continue;
        const checked = positive(value, label);
        if (typeof checked !== 'number')
            return checked;
    }
    const missing = ['pressureAtm', 'volumeL', 'moles', 'temperatureK'].find((key) => input[key] === undefined);
    if (missing === 'pressureAtm') {
        const moles = input.moles;
        const temperatureK = input.temperatureK;
        const volumeL = input.volumeL;
        const pressureAtm = (moles * R_ATM_L * temperatureK) / volumeL;
        return { ok: true, solved: 'pressure', ...withConversions(pressureAtm, volumeL, moles, temperatureK) };
    }
    if (missing === 'volumeL') {
        const moles = input.moles;
        const temperatureK = input.temperatureK;
        const pressureAtm = input.pressureAtm;
        const volumeL = (moles * R_ATM_L * temperatureK) / pressureAtm;
        return { ok: true, solved: 'volume', ...withConversions(pressureAtm, volumeL, moles, temperatureK) };
    }
    if (missing === 'moles') {
        const pressureAtm = input.pressureAtm;
        const volumeL = input.volumeL;
        const temperatureK = input.temperatureK;
        const moles = (pressureAtm * volumeL) / (R_ATM_L * temperatureK);
        return { ok: true, solved: 'moles', ...withConversions(pressureAtm, volumeL, moles, temperatureK) };
    }
    const pressureAtm = input.pressureAtm;
    const volumeL = input.volumeL;
    const moles = input.moles;
    const temperatureK = (pressureAtm * volumeL) / (moles * R_ATM_L);
    return { ok: true, solved: 'temperature', ...withConversions(pressureAtm, volumeL, moles, temperatureK) };
}
/** Molar volume V/n = RT/P (defaults: 0 °C, 1 atm → 22.414 L/mol). */
export function molarVolume(temperatureK = 273.15, pressureAtm = 1) {
    const t = positive(temperatureK, 'temperatureK');
    if (typeof t !== 'number')
        return t;
    const p = positive(pressureAtm, 'pressureAtm');
    if (typeof p !== 'number')
        return p;
    const perMol = (R_ATM_L * t) / p;
    return {
        ok: true,
        temperatureK: t,
        temperatureC: t - KELVIN_OFFSET,
        pressureAtm: p,
        pressureKpa: p * ATM_TO_KPA,
        molarVolumeLPerMol: perMol,
        molarVolumeMlPerMol: perMol * 1000,
        r: R_ATM_L,
    };
}
/** Gas density from the ideal gas law: ρ = P·M / (R·T). */
export function gasDensity(pressureAtm, molarMass, temperatureK) {
    const p = positive(pressureAtm, 'pressureAtm');
    if (typeof p !== 'number')
        return p;
    const m = positive(molarMass, 'molarMass');
    if (typeof m !== 'number')
        return m;
    const t = positive(temperatureK, 'temperatureK');
    if (typeof t !== 'number')
        return t;
    const density = (p * m) / (R_ATM_L * t);
    return {
        ok: true,
        pressureAtm: p,
        temperatureK: t,
        temperatureC: t - KELVIN_OFFSET,
        molarMass: m,
        densityGPerL: density,
        densityGPerMl: density / 1000,
        r: R_ATM_L,
    };
}
/** Molar mass back-calculated from a measured gas density: M = ρ·R·T / P. */
export function molarMassFromDensity(densityGPerL, pressureAtm, temperatureK) {
    const d = positive(densityGPerL, 'densityGPerL');
    if (typeof d !== 'number')
        return d;
    const p = positive(pressureAtm, 'pressureAtm');
    if (typeof p !== 'number')
        return p;
    const t = positive(temperatureK, 'temperatureK');
    if (typeof t !== 'number')
        return t;
    return {
        ok: true,
        pressureAtm: p,
        temperatureK: t,
        temperatureC: t - KELVIN_OFFSET,
        densityGPerL: d,
        molarMass: (d * R_ATM_L * t) / p,
        r: R_ATM_L,
    };
}
//# sourceMappingURL=gaslaw.js.map