"""Build the nuclear reaction networks used by the nova hydrogen-burning study.

This script is the only step that needs an external nuclear-data source.  It
reads the JINA ReacLib snapshot and the atomic mass evaluation that ship with
``pynucastro``, selects the nuclides belonging to each network case of the
proposal (Table "Planned nuclear network cases"), and writes each network out
in two forms:

``data/networks/<case>.json``
    A complete, self-describing archive: every nuclide with its charge, mass
    number and mass excess, and every reaction with its ReacLib seven-parameter
    coefficient sets, Q value, source label and reverse flag.  This is what the
    experiment reads, so the calculations can be repeated later with nothing
    but NucNetPy, NumPy and SciPy installed.

``data/networks/<case>.xml``
    The same network in NucNetPy's own XML format, so that it can also be
    inspected with the ``nucnetpy`` command-line tool.

Run it once::

    python src/build_networks.py

Nuclide selection
-----------------
For every network case we keep the nuclides that ReacLib knows about with
``Z <= Z_max`` whose neutron excess satisfies ``N - Z <= 3 + Z/4``.  Nova
hydrogen burning runs along the proton-rich side of the valley of stability, so
nuclides beyond that band never acquire an abundance; leaving them out keeps
the integration affordable without changing the result, and the cut is wide
enough to retain every stable isotope of every element in the starting
composition.  Neutrons, protons, deuterons, tritons, :sup:`3`\\ He and
:sup:`4`\\ He are always kept.  A few nuclides beyond the proton drip line have
no measured mass excess and are dropped as well.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NETWORK_DIR = ROOT / "data" / "networks"

#: Network cases of the proposal: label -> (maximum charge, description).
NETWORK_CASES = {
    "nova_z10": (10, "Small network, Z <= 10: CNO-limited calculation"),
    "nova_z20": (20, "Intermediate network, Z <= 20: extension beyond CNO"),
    "nova_z30": (30, "Large network, Z <= 30: extended nova network"),
}

#: Neutron-richness cut, as a limit on ``N - Z``.
#:
#: The limit has to widen with charge, because the valley of stability does:
#: 40Ca sits at ``N = Z`` but 56Fe sits four neutrons beyond it and 58Fe six.
#: A flat cut would throw away the most abundant iron isotope in the starting
#: composition.  ``N - Z <= 3 + Z/4`` keeps every stable isotope of every
#: element in the solar composition, with room to spare on the neutron-rich
#: side, while still leaving out the far side of the chart, which a
#: proton-rich nova envelope can never populate.
def max_neutron_excess(z: int) -> float:
    return 3.0 + z / 4.0

#: Light nuclides that are always part of the network.
ALWAYS_KEEP = {"n", "p", "d", "t", "he3", "he4"}

#: Particle-unbound nuclides, with what they break up into.
#:
#: ReacLib produces these three but gives them no decay channel, because they
#: fall apart in less than 10^-16 s and no network is meant to carry them as
#: abundances.  Left in a network they become dead ends that quietly swallow
#: material.  Every reaction that makes one is therefore rewritten to make its
#: break-up products instead, which is how such decays are normally written
#: (8B beta-decays to 2 alpha, not to 8Be).  None of the three ever appears as
#: a reactant in the library, so the substitution is complete.
UNBOUND_PRODUCTS = {
    "li5": ["he4", "p"],
    "be8": ["he4", "he4"],
    "b9": ["he4", "he4", "p"],
}

#: A one-body rate faster than this (per second, evaluated at the coldest
#: temperature the calculations reach) does not describe a nuclide that a
#: network can carry: it means prompt particle emission.  ReacLib does give
#: such nuclides a decay -- 18Na -> p + 17Ne runs at 5e20 s^-1 -- and keeping
#: them makes the system stiff by twenty orders of magnitude for no physics.
#: They are eliminated the same way as the dead ends above, by substituting
#: their break-up products wherever they are produced.
PROMPT_DECAY_RATE = 1.0e12

#: Temperature at which the prompt-decay test is applied.  It is the lower end
#: of the range the calculations use, and it is cold enough that photo-
#: disintegration cannot make a bound nuclide look unbound.
TEST_T9 = 0.01


def _keep_nucleus(z: int, a: int, raw: str) -> bool:
    if raw in ALWAYS_KEEP:
        return True
    return (a - 2 * z) <= max_neutron_excess(z)


def _reaclib_rate(sets: list[list[float]], t9: float) -> float:
    """Evaluate summed ReacLib fits, the same expression NucNetPy uses."""
    total = 0.0
    for a in sets:
        t13 = t9 ** (1.0 / 3.0)
        exponent = (a[0] + a[1] / t9 + a[2] / t13 + a[3] * t13
                    + a[4] * t9 + a[5] * t9 ** (5.0 / 3.0) + a[6] * math.log(t9))
        total += math.exp(exponent) if exponent < 700.0 else float("inf")
    return total


def find_prompt_emitters(reactions: list[dict]) -> dict:
    """Nuclides that fall apart promptly, mapped to what they fall apart into."""
    fastest: dict[str, tuple] = {}
    for entry in reactions:
        if len(entry["reactants"]) != 1:
            continue
        rate = _reaclib_rate(entry["sets"], TEST_T9)
        if rate <= PROMPT_DECAY_RATE:
            continue
        parent = entry["reactants"][0]
        if parent not in fastest or rate > fastest[parent][0]:
            fastest[parent] = (rate, list(entry["products"]))
    return {parent: products for parent, (_, products) in fastest.items()}


def substitute_unbound(reactions: list[dict], unbound: dict) -> list[dict]:
    """Rewrite every reaction so that no unbound nuclide appears in it.

    Products are replaced by their break-up products, repeatedly, in case one
    break-up produces another unbound nuclide.  Reactions that have an unbound
    nuclide as a reactant are dropped, because that nuclide no longer exists;
    so are reactions whose products end up identical to their reactants, which
    is what a capture onto an unbound nuclide becomes once the prompt
    re-emission is folded in.
    """
    out = []
    for entry in reactions:
        if any(name in unbound for name in entry["reactants"]):
            continue

        products = list(entry["products"])
        changed = False
        while any(name in unbound for name in products):
            expanded = []
            for name in products:
                expanded.extend(unbound.get(name, [name]))
            products = expanded
            changed = True

        if sorted(products) == sorted(entry["reactants"]):
            continue

        entry = dict(entry)
        entry["products"] = products
        entry["substituted"] = changed
        out.append(entry)
    return out


def build_case(library, z_max: int) -> dict:
    """Return a network archive dictionary for one ``Z <= z_max`` case."""
    from pynucastro.nucdata import Nucleus

    nuclides: dict[str, Nucleus] = {}
    reactions: list[dict] = []
    seen: set[tuple] = set()

    for rate in library.get_rates():
        parts = list(rate.reactants) + list(rate.products)
        if any(nuc.Z > z_max for nuc in parts):
            continue
        if not all(_keep_nucleus(nuc.Z, nuc.A, nuc.raw) for nuc in parts):
            continue
        # A handful of nuclides beyond the proton drip line have no measured
        # mass excess.  They cannot be populated at nova temperatures, and
        # without a mass they would break the energy-generation bookkeeping, so
        # they and their reactions are dropped.
        if any(nuc.dm is None for nuc in parts):
            continue

        # ReacLib occasionally carries two evaluations of the same reaction.
        # Keep the first one encountered so that the network stays a single
        # consistent rate per reaction rather than a double-counted sum.
        key = (
            tuple(sorted(nuc.raw for nuc in rate.reactants)),
            tuple(sorted(nuc.raw for nuc in rate.products)),
        )
        if key in seen:
            continue
        seen.add(key)

        for nuc in parts:
            nuclides.setdefault(nuc.raw, nuc)

        reactions.append(
            {
                "reactants": [nuc.raw for nuc in rate.reactants],
                "products": [nuc.raw for nuc in rate.products],
                "q_value": float(rate.Q),
                "label": str(rate.label).strip(),
                "chapter": int(rate.chapter) if str(rate.chapter).isdigit() else None,
                "weak": bool(rate.weak),
                "reverse": str(rate.labelprops).strip().endswith("v"),
                "sets": [[float(c) for c in s.a] for s in rate.sets],
            }
        )

    # Remove the nuclides that no network can carry: the dead ends the library
    # gives no decay at all, and the prompt particle emitters it gives an
    # absurdly fast one.
    unbound = dict(UNBOUND_PRODUCTS)
    unbound.update(find_prompt_emitters(reactions))
    reactions = substitute_unbound(reactions, unbound)
    for name in unbound:
        nuclides.pop(name, None)

    # Rebuild the Q value of every rewritten reaction from the mass excesses of
    # the products actually written down; the library value stopped at the
    # unbound nuclide.
    for entry in reactions:
        if entry.pop("substituted", False):
            entry["q_value"] = (
                sum(Nucleus(name).dm for name in entry["reactants"])
                - sum(Nucleus(name).dm for name in entry["products"])
            )

    species = [
        {
            "name": nuc.raw,
            "z": int(nuc.Z),
            "a": int(nuc.A),
            "mass_excess": float(nuc.dm),
            "spin": (float(nuc.spin_states) - 1.0) / 2.0 if nuc.spin_states else None,
        }
        for nuc in sorted(nuclides.values(), key=lambda n: (n.Z, n.A))
    ]

    reactions.sort(key=lambda r: (sorted(r["reactants"]), sorted(r["products"])))
    return {
        "z_max": z_max,
        "eliminated_nuclides": {k: v for k, v in sorted(unbound.items())},
        "species": species,
        "reactions": reactions,
    }


def main() -> None:
    import pynucastro as pyna

    NETWORK_DIR.mkdir(parents=True, exist_ok=True)
    library = pyna.ReacLibLibrary()

    # Imported lazily so that the rest of the experiment never needs NucNetPy's
    # XML writer unless the networks are rebuilt.
    from network_io import archive_to_network
    from nucnetpy import write_xml

    for case, (z_max, description) in NETWORK_CASES.items():
        archive = build_case(library, z_max)
        archive["description"] = description
        archive["source"] = "JINA ReacLib snapshot 20180319default2 (via pynucastro)"
        archive["mass_source"] = "AME atomic mass evaluation (via pynucastro)"

        json_path = NETWORK_DIR / f"{case}.json"
        json_path.write_text(json.dumps(archive, indent=1))

        write_xml(archive_to_network(archive), NETWORK_DIR / f"{case}.xml")

        print(
            f"{case:10s} Z<={z_max:2d}  "
            f"{len(archive['species']):4d} nuclides  "
            f"{len(archive['reactions']):5d} reactions  -> {json_path.name}"
        )


if __name__ == "__main__":
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
