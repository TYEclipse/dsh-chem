/**
 * Periodic table data for dsh-chem.
 *
 * Atomic masses are the IUPAC standard atomic weights (abridged, conventional
 * values, widely used in textbooks). Elements with no stable isotopes use the
 * mass number of the most stable known isotope (marked with `massNumber`),
 * matching the bracketed-value convention of IUPAC tables.
 *
 * @module dsh-chem/periodic
 */

export interface ElementData {
  /** Element symbol, case-sensitive ("Co" is cobalt, "CO" is carbon monoxide). */
  symbol: string
  name: string
  /** Atomic number Z (1–118). */
  atomicNumber: number
  /** Standard atomic weight in g/mol. */
  atomicMass: number
  /** True when the mass is the most-stable-isotope mass number (radioelement). */
  massNumber: boolean
}

const el = (symbol: string, name: string, atomicNumber: number, atomicMass: number, massNumber = false): ElementData =>
  ({ symbol, name, atomicNumber, atomicMass, massNumber })

/** All 118 elements, ordered by atomic number. */
export const ELEMENTS: readonly ElementData[] = [
  el('H', 'Hydrogen', 1, 1.008),
  el('He', 'Helium', 2, 4.0026),
  el('Li', 'Lithium', 3, 6.94),
  el('Be', 'Beryllium', 4, 9.0122),
  el('B', 'Boron', 5, 10.81),
  el('C', 'Carbon', 6, 12.011),
  el('N', 'Nitrogen', 7, 14.007),
  el('O', 'Oxygen', 8, 15.999),
  el('F', 'Fluorine', 9, 18.998),
  el('Ne', 'Neon', 10, 20.18),
  el('Na', 'Sodium', 11, 22.99),
  el('Mg', 'Magnesium', 12, 24.305),
  el('Al', 'Aluminium', 13, 26.982),
  el('Si', 'Silicon', 14, 28.085),
  el('P', 'Phosphorus', 15, 30.974),
  el('S', 'Sulfur', 16, 32.06),
  el('Cl', 'Chlorine', 17, 35.45),
  el('Ar', 'Argon', 18, 39.948),
  el('K', 'Potassium', 19, 39.098),
  el('Ca', 'Calcium', 20, 40.078),
  el('Sc', 'Scandium', 21, 44.956),
  el('Ti', 'Titanium', 22, 47.867),
  el('V', 'Vanadium', 23, 50.942),
  el('Cr', 'Chromium', 24, 51.996),
  el('Mn', 'Manganese', 25, 54.938),
  el('Fe', 'Iron', 26, 55.845),
  el('Co', 'Cobalt', 27, 58.933),
  el('Ni', 'Nickel', 28, 58.693),
  el('Cu', 'Copper', 29, 63.546),
  el('Zn', 'Zinc', 30, 65.38),
  el('Ga', 'Gallium', 31, 69.723),
  el('Ge', 'Germanium', 32, 72.63),
  el('As', 'Arsenic', 33, 74.922),
  el('Se', 'Selenium', 34, 78.971),
  el('Br', 'Bromine', 35, 79.904),
  el('Kr', 'Krypton', 36, 83.798),
  el('Rb', 'Rubidium', 37, 85.468),
  el('Sr', 'Strontium', 38, 87.62),
  el('Y', 'Yttrium', 39, 88.906),
  el('Zr', 'Zirconium', 40, 91.224),
  el('Nb', 'Niobium', 41, 92.906),
  el('Mo', 'Molybdenum', 42, 95.95),
  el('Tc', 'Technetium', 43, 98, true),
  el('Ru', 'Ruthenium', 44, 101.07),
  el('Rh', 'Rhodium', 45, 102.91),
  el('Pd', 'Palladium', 46, 106.42),
  el('Ag', 'Silver', 47, 107.87),
  el('Cd', 'Cadmium', 48, 112.41),
  el('In', 'Indium', 49, 114.82),
  el('Sn', 'Tin', 50, 118.71),
  el('Sb', 'Antimony', 51, 121.76),
  el('Te', 'Tellurium', 52, 127.6),
  el('I', 'Iodine', 53, 126.9),
  el('Xe', 'Xenon', 54, 131.29),
  el('Cs', 'Caesium', 55, 132.91),
  el('Ba', 'Barium', 56, 137.33),
  el('La', 'Lanthanum', 57, 138.91),
  el('Ce', 'Cerium', 58, 140.12),
  el('Pr', 'Praseodymium', 59, 140.91),
  el('Nd', 'Neodymium', 60, 144.24),
  el('Pm', 'Promethium', 61, 145, true),
  el('Sm', 'Samarium', 62, 150.36),
  el('Eu', 'Europium', 63, 151.96),
  el('Gd', 'Gadolinium', 64, 157.25),
  el('Tb', 'Terbium', 65, 158.93),
  el('Dy', 'Dysprosium', 66, 162.5),
  el('Ho', 'Holmium', 67, 164.93),
  el('Er', 'Erbium', 68, 167.26),
  el('Tm', 'Thulium', 69, 168.93),
  el('Yb', 'Ytterbium', 70, 173.05),
  el('Lu', 'Lutetium', 71, 174.97),
  el('Hf', 'Hafnium', 72, 178.49),
  el('Ta', 'Tantalum', 73, 180.95),
  el('W', 'Tungsten', 74, 183.84),
  el('Re', 'Rhenium', 75, 186.21),
  el('Os', 'Osmium', 76, 190.23),
  el('Ir', 'Iridium', 77, 192.22),
  el('Pt', 'Platinum', 78, 195.08),
  el('Au', 'Gold', 79, 196.97),
  el('Hg', 'Mercury', 80, 200.59),
  el('Tl', 'Thallium', 81, 204.38),
  el('Pb', 'Lead', 82, 207.2),
  el('Bi', 'Bismuth', 83, 208.98),
  el('Po', 'Polonium', 84, 209, true),
  el('At', 'Astatine', 85, 210, true),
  el('Rn', 'Radon', 86, 222, true),
  el('Fr', 'Francium', 87, 223, true),
  el('Ra', 'Radium', 88, 226, true),
  el('Ac', 'Actinium', 89, 227, true),
  el('Th', 'Thorium', 90, 232.04),
  el('Pa', 'Protactinium', 91, 231.04),
  el('U', 'Uranium', 92, 238.03),
  el('Np', 'Neptunium', 93, 237, true),
  el('Pu', 'Plutonium', 94, 244, true),
  el('Am', 'Americium', 95, 243, true),
  el('Cm', 'Curium', 96, 247, true),
  el('Bk', 'Berkelium', 97, 247, true),
  el('Cf', 'Californium', 98, 251, true),
  el('Es', 'Einsteinium', 99, 252, true),
  el('Fm', 'Fermium', 100, 257, true),
  el('Md', 'Mendelevium', 101, 258, true),
  el('No', 'Nobelium', 102, 259, true),
  el('Lr', 'Lawrencium', 103, 266, true),
  el('Rf', 'Rutherfordium', 104, 267, true),
  el('Db', 'Dubnium', 105, 268, true),
  el('Sg', 'Seaborgium', 106, 269, true),
  el('Bh', 'Bohrium', 107, 270, true),
  el('Hs', 'Hassium', 108, 277, true),
  el('Mt', 'Meitnerium', 109, 278, true),
  el('Ds', 'Darmstadtium', 110, 281, true),
  el('Rg', 'Roentgenium', 111, 282, true),
  el('Cn', 'Copernicium', 112, 285, true),
  el('Nh', 'Nihonium', 113, 286, true),
  el('Fl', 'Flerovium', 114, 289, true),
  el('Mc', 'Moscovium', 115, 290, true),
  el('Lv', 'Livermorium', 116, 293, true),
  el('Ts', 'Tennessine', 117, 294, true),
  el('Og', 'Oganesson', 118, 294, true),
]

/** Lookup by exact case-sensitive symbol. */
const BY_SYMBOL: ReadonlyMap<string, ElementData> = new Map(ELEMENTS.map((e) => [e.symbol, e]))

/** Lookup by case-insensitive symbol (fallback for lenient matching). */
const BY_SYMBOL_LOWER: ReadonlyMap<string, ElementData> = new Map(ELEMENTS.map((e) => [e.symbol.toLowerCase(), e]))

export function elementBySymbol(symbol: string): ElementData | undefined {
  return BY_SYMBOL.get(symbol)
}

/**
 * Lenient lookup: exact symbol first, then lowercase-insensitive match
 * ("co" -> Co cobalt). This catches the most common LLM mistake (lowercasing
 * symbols) while still treating "CO" as carbon + oxygen.
 */
export function elementBySymbolLenient(symbol: string): ElementData | undefined {
  const exact = BY_SYMBOL.get(symbol)
  if (exact) return exact
  return BY_SYMBOL_LOWER.get(symbol.toLowerCase())
}
