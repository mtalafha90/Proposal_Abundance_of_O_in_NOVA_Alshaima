"""Turn a stored network archive into NucNetPy objects.

The archives in ``data/networks`` hold the nuclear data in a plain JSON form
(see ``build_networks.py``).  This module is the bridge from that data to the
NucNetPy classes that actually do the physics: :class:`nucnetpy.Species`,
:class:`nucnetpy.Reaction` and :class:`nucnetpy.Network`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from nucnetpy import Network, RateFit, Reaction, Species
from nucnetpy.network_limiter import limit_network, select_species

ROOT = Path(__file__).resolve().parent.parent
NETWORK_DIR = ROOT / "data" / "networks"


def load_archive(case: str) -> dict:
    """Read the JSON archive for a network case such as ``"nova_z10"``."""
    return json.loads((NETWORK_DIR / f"{case}.json").read_text())


def archive_to_network(archive: dict) -> Network:
    """Build a :class:`nucnetpy.Network` from a stored archive."""
    net = Network()
    for entry in archive["species"]:
        net.add_species(
            Species.parse(
                entry["name"],
                mass_excess=entry["mass_excess"],
                spin=entry.get("spin"),
                source=archive.get("mass_source", ""),
            )
        )

    for entry in archive["reactions"]:
        # Each ReacLib "set" is one additive term of the same rate (a resonant
        # and a non-resonant contribution, say), which is exactly how NucNetPy
        # treats a list of rate fits: it sums them.
        fits = [RateFit(coeffs, label=entry.get("label", "")) for coeffs in entry["sets"]]
        net.reactions.add(
            Reaction.from_names(
                entry["reactants"],
                entry["products"],
                rate_fits=fits,
                q_value=entry.get("q_value", 0.0),
                label=entry.get("label", ""),
                source=entry.get("label", ""),
                metadata={
                    "weak": str(bool(entry.get("weak"))),
                    "reverse": str(bool(entry.get("reverse"))),
                },
            )
        )
    return net


def load_network(case: str) -> Network:
    """Load a network case by name."""
    return archive_to_network(load_archive(case))


def restrict_to_charge(net: Network, z_max: int) -> Network:
    """Cut a network down to nuclides with ``Z <= z_max`` using NucNetPy."""
    return limit_network(net, select_species(net, zmax=z_max))


def find_reaction(net: Network, reactants: Iterable[str], products: Iterable[str]):
    """Return the reaction with the given reactants and products, or ``None``.

    Reactions are looked up by what they do rather than by label, because the
    NucNetPy XML round trip does not preserve labels.
    """
    want = (
        tuple(sorted(Species.parse(r).name if r != "gamma" else "gamma" for r in reactants)),
        tuple(sorted(Species.parse(p).name if p != "gamma" else "gamma" for p in products)),
    )
    for reaction in net.reactions.reactions:
        key = (
            tuple(sorted(p.species for p in reaction.reactants for _ in range(p.count))),
            tuple(sorted(p.species for p in reaction.products for _ in range(p.count))),
        )
        if key == want:
            return reaction
    return None
