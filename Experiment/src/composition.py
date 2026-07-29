"""Solar initial composition for the nova hydrogen-burning calculations.

The proposal starts every calculation from a solar-like composition taken from
Bergemann, Lodders & Palme (2025).  This module encodes that composition as
elemental mass fractions together with the solar isotopic splits, and converts
it into the molar abundances ``Y_i = X_i / A_i`` that NucNetPy evolves.

The two elements that define the diagnostic ratio are pinned to the values
quoted in the proposal:

    X(14N) = 7.97e-4      X(15N) = 3.14e-6

which follow from a nitrogen mass fraction of 8.00e-4 and the solar isotopic
ratio 15N/14N = 1/272, and give an initial diagnostic ratio

    R_15/14 = 3.9e-3.

Hydrogen carries whatever mass fraction is left over, so the composition sums
to one by construction.
"""

from __future__ import annotations

from typing import Dict

from nucnetpy import Species

#: Elemental mass fractions of the solar composition (Bergemann, Lodders &
#: Palme 2025).  Hydrogen is not listed: it takes up the remainder.
ELEMENT_MASS_FRACTIONS: Dict[str, float] = {
    "he": 2.757e-1,
    "c": 2.30e-3,
    "n": 8.00e-4,
    "o": 5.79e-3,
    "ne": 1.55e-3,
    "na": 3.40e-5,
    "mg": 6.50e-4,
    "al": 5.60e-5,
    "si": 6.70e-4,
    "s": 3.40e-4,
    "ar": 8.70e-5,
    "ca": 6.40e-5,
    "fe": 1.30e-3,
}

#: Solar isotopic composition, given as number fractions within each element.
#: Values are the standard solar-system isotopic abundances; each list is
#: normalised internally, so small rounding in the entries does not matter.
ISOTOPE_NUMBER_FRACTIONS: Dict[str, Dict[int, float]] = {
    # Deuterium is set to zero rather than to its protosolar value.  Material
    # accreted onto the white dwarf has already been through the companion
    # star, where deuterium is destroyed long before the main sequence.  Left
    # in, it burns away in the first microsecond of the calculation and puts a
    # spurious spike on the energy-generation rate; it does not touch the CNO
    # isotopes either way.
    "h": {1: 1.0},
    "he": {3: 1.66e-4, 4: 1.0 - 1.66e-4},
    "c": {12: 0.98890, 13: 0.01110},
    "n": {14: 1.0 / (1.0 + 1.0 / 272.0), 15: (1.0 / 272.0) / (1.0 + 1.0 / 272.0)},
    "o": {16: 0.997620, 17: 0.000380, 18: 0.002000},
    "ne": {20: 0.92940, 21: 0.00222, 22: 0.06838},
    "na": {23: 1.0},
    "mg": {24: 0.78990, 25: 0.10000, 26: 0.11010},
    "al": {27: 1.0},
    "si": {28: 0.92230, 29: 0.04670, 30: 0.03100},
    "s": {32: 0.94930, 33: 0.00760, 34: 0.04290, 36: 0.00020},
    "ar": {36: 0.84600, 38: 0.15400},
    "ca": {40: 0.96940, 42: 0.00647, 43: 0.00135, 44: 0.02090},
    "fe": {54: 0.05845, 56: 0.91754, 57: 0.02119, 58: 0.00282},
}


def solar_mass_fractions() -> Dict[str, float]:
    """Return the solar composition as ``{species name: mass fraction}``."""
    heavy: Dict[str, float] = {}
    for element, x_element in ELEMENT_MASS_FRACTIONS.items():
        splits = ISOTOPE_NUMBER_FRACTIONS[element]
        mean_a = sum(f * a for a, f in splits.items()) / sum(splits.values())
        for a, number_fraction in splits.items():
            weight = number_fraction * a / (mean_a * sum(splits.values()))
            heavy[f"{element}{a}"] = x_element * weight

    x_hydrogen = 1.0 - sum(heavy.values())
    splits = ISOTOPE_NUMBER_FRACTIONS["h"]
    mean_a = sum(f * a for a, f in splits.items())
    out = {f"h{a}": x_hydrogen * f * a / mean_a for a, f in splits.items()}
    out.update(heavy)
    return out


def solar_abundances(species_names=None) -> Dict[str, float]:
    """Return the solar composition as molar abundances ``Y = X / A``.

    If ``species_names`` is given, species outside that set are dropped and the
    composition is renormalised so that ``sum(A_i Y_i) = 1`` is preserved.  A
    small network cannot hold the whole solar composition -- a ``Z <= 10``
    network has nowhere to put iron -- so this is expected rather than an
    error, but how much was dropped should always be reported alongside the
    result.  Use :func:`truncation_report` for that.  Note that renormalising
    scales every surviving abundance by the same factor, so it cannot change
    the diagnostic ratio, which is a ratio of abundances.
    """
    mass_fractions = solar_mass_fractions()
    if species_names is not None:
        allowed = set(species_names)
        mass_fractions = {k: v for k, v in mass_fractions.items() if k in allowed}
        total = sum(mass_fractions.values())
        mass_fractions = {k: v / total for k, v in mass_fractions.items()}
    return {name: x / Species.parse(name).a for name, x in mass_fractions.items()}


def truncation_report(species_names) -> Dict[str, object]:
    """What the starting composition loses when restricted to a network."""
    allowed = set(species_names)
    dropped = {k: v for k, v in solar_mass_fractions().items() if k not in allowed}
    return {
        "dropped_species": sorted(dropped),
        "dropped_mass_fraction": sum(dropped.values()),
        "renormalisation_factor": 1.0 / (1.0 - sum(dropped.values())),
    }
