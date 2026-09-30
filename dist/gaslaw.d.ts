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
export declare const R_JOULES = 8.31446261815324;
/** 1 atm in kPa (exact by definition of the standard atmosphere). */
export declare const ATM_TO_KPA = 101.325;
/** 1 atm in mmHg (conventional). */
export declare const ATM_TO_MMHG = 760;
/** Molar gas constant in L·atm/(mol·K). */
export declare const R_ATM_L: number;
/** Kelvin offset from degrees Celsius. */
export declare const KELVIN_OFFSET = 273.15;
export interface GasQuantities {
    pressureAtm?: number | undefined;
    volumeL?: number | undefined;
    moles?: number | undefined;
    temperatureK?: number | undefined;
}
export type GasFailure = {
    ok: false;
    error: string;
};
export interface IdealGasSolution {
    ok: true;
    /** Which of the four quantities was computed from the other three. */
    solved: 'pressure' | 'volume' | 'moles' | 'temperature';
    pressureAtm: number;
    pressureKpa: number;
    pressureMmhg: number;
    volumeL: number;
    volumeMl: number;
    moles: number;
    temperatureK: number;
    temperatureC: number;
    r: number;
}
export interface MolarVolumeSolution {
    ok: true;
    temperatureK: number;
    temperatureC: number;
    pressureAtm: number;
    pressureKpa: number;
    molarVolumeLPerMol: number;
    molarVolumeMlPerMol: number;
    r: number;
}
export interface GasDensitySolution {
    ok: true;
    pressureAtm: number;
    temperatureK: number;
    temperatureC: number;
    molarMass: number;
    densityGPerL: number;
    densityGPerMl: number;
    r: number;
}
export interface GasMolarMassSolution {
    ok: true;
    pressureAtm: number;
    temperatureK: number;
    temperatureC: number;
    densityGPerL: number;
    molarMass: number;
    r: number;
}
/** Solve PV = nRT for whichever of the four quantities was left out. */
export declare function idealGas(input: GasQuantities): IdealGasSolution | GasFailure;
/** Molar volume V/n = RT/P (defaults: 0 °C, 1 atm → 22.414 L/mol). */
export declare function molarVolume(temperatureK?: number, pressureAtm?: number): MolarVolumeSolution | GasFailure;
/** Gas density from the ideal gas law: ρ = P·M / (R·T). */
export declare function gasDensity(pressureAtm: number, molarMass: number, temperatureK: number): GasDensitySolution | GasFailure;
/** Molar mass back-calculated from a measured gas density: M = ρ·R·T / P. */
export declare function molarMassFromDensity(densityGPerL: number, pressureAtm: number, temperatureK: number): GasMolarMassSolution | GasFailure;
//# sourceMappingURL=gaslaw.d.ts.map