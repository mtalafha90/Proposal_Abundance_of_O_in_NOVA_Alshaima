"""Check the networks before they are used, and record the checks.

Run it with::

    python src/validate_networks.py

It writes ``results/network_validation.json`` and prints a summary.  Six things
are checked for each network case:

1.  **Baryon number.**  Every reaction must have ``sum(A)`` equal on both sides.
    NucNetPy's own ``conserves_a_z`` also demands charge balance, which weak
    decays fail only because ReacLib leaves the leptons out of the record; so
    charge is reported separately, and the only acceptable charge changes are
    the ``+-1`` of a beta decay or an electron capture.

2.  **No dead ends.**  Every nuclide that a reaction can produce must also have
    a way out, otherwise it silently accumulates material.

3.  **Usable rates.**  Over the whole temperature range the calculations visit,
    every rate must be finite, and no one-body rate may exceed
    ``1e12 per second``: a nuclide that decays that fast is not one a network
    can carry, and it makes the system stiff for no physics.

4.  **The XML export.**  The NucNetPy XML written beside each archive must
    describe the same network: the same nuclide and reaction counts, the same
    mass excesses, and the same rates.

5.  **Known rates.**  The rates of the reactions that control the diagnostic
    ratio are compared against their measured half-lives, which is an
    end-to-end test of the fitting function, the coefficients and the units.

6.  **The composition.**  The initial mass fractions must sum to one, and the
    two isotopes that define the diagnostic ratio must match the values quoted
    in the proposal.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import composition
import thermodynamics
from network_io import find_reaction, load_network

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

CASES = ["nova_z10", "nova_z20", "nova_z30"]

#: Beta-decay half-lives in seconds, for the end-to-end rate check.
HALF_LIVES = {
    ("n13", "c13"): 597.9,
    ("o14", "n14"): 70.62,
    ("o15", "n15"): 122.24,
    ("f17", "o17"): 64.37,
    ("f18", "o18"): 6586.0,
}


def check_conservation(net) -> dict:
    baryon_violations = []
    charge_changes = {}
    for reaction in net.reactions.reactions:
        delta_a = delta_z = 0
        for name, nu in reaction.stoichiometry().items():
            species = net.species.get(name)
            if species is None:
                continue
            delta_a += nu * species.a
            delta_z += nu * species.z
        if delta_a != 0:
            baryon_violations.append(reaction.string)
        if delta_z != 0:
            charge_changes[delta_z] = charge_changes.get(delta_z, 0) + 1
    return {
        "baryon_violations": baryon_violations,
        "charge_changes": {str(k): v for k, v in sorted(charge_changes.items())},
        "charge_changes_are_all_beta": set(charge_changes) <= {-1, 1},
    }


def check_dead_ends(net) -> list:
    produced, consumed = set(), set()
    for reaction in net.reactions.reactions:
        for part in reaction.products:
            produced.add(part.species)
        for part in reaction.reactants:
            consumed.add(part.species)
    return sorted(produced - consumed)


def check_rate_sanity(net) -> dict:
    """No infinite rates, and no prompt one-body decays, anywhere in range."""
    temperatures = [thermodynamics.T9_FLOOR, 0.05, 0.1, 0.2, 0.45, 1.0]
    non_finite, prompt = [], []
    for reaction in net.reactions.reactions:
        for t9 in temperatures:
            rate = reaction.rate(t9)
            if not math.isfinite(rate):
                non_finite.append([reaction.string, t9])
                break
            if reaction.reactant_order == 1 and rate > 1.0e12:
                prompt.append([reaction.string, t9, rate])
                break
    return {"non_finite": non_finite, "prompt_one_body_decays": prompt}


def check_rates(net) -> dict:
    out = {}
    for (parent, daughter), half_life in HALF_LIVES.items():
        reaction = find_reaction(net, [parent], [daughter])
        if reaction is None:
            continue
        # A beta decay is temperature independent, so any T9 will do.
        rate = reaction.rate(0.1)
        measured = math.log(2.0) / half_life
        out[f"{parent}->{daughter}"] = {
            "rate_from_network": rate,
            "rate_from_half_life": measured,
            "relative_difference": (rate - measured) / measured,
        }
    return out


def check_xml_export(case: str, net) -> dict:
    """The NucNetPy XML export must describe the same network as the archive."""
    from nucnetpy import read_xml
    from network_io import NETWORK_DIR

    other = read_xml(NETWORK_DIR / f"{case}.xml")
    rates_json = sorted(r.rate(0.2) for r in net.reactions.reactions)
    rates_xml = sorted(r.rate(0.2) for r in other.reactions.reactions)
    masses_json = sorted(s.mass_excess for s in net.species.values())
    masses_xml = sorted(s.mass_excess for s in other.species.values())
    same_size = (len(rates_json) == len(rates_xml)
                 and len(masses_json) == len(masses_xml))
    return {
        "same_counts": same_size,
        "largest_rate_difference": (
            max(abs(a - b) for a, b in zip(rates_json, rates_xml)) if same_size else None
        ),
        "largest_mass_difference": (
            max(abs(a - b) for a, b in zip(masses_json, masses_xml)) if same_size else None
        ),
    }


def check_composition() -> dict:
    mass_fractions = composition.solar_mass_fractions()
    abundances = composition.solar_abundances()
    ratio = (abundances["n15"]) / (abundances["n14"])
    return {
        "sum_of_mass_fractions": sum(mass_fractions.values()),
        "x_n14": mass_fractions["n14"],
        "x_n15": mass_fractions["n15"],
        "x_n14_in_proposal": 7.97e-4,
        "x_n15_in_proposal": 3.14e-6,
        "r_initial_molar": ratio,
        "r_initial_by_mass_fraction": mass_fractions["n15"] / mass_fractions["n14"],
    }


def main() -> None:
    report = {"composition": check_composition(), "networks": {}}
    for case in CASES:
        net = load_network(case)
        report["networks"][case] = {
            "species": len(net.species),
            "reactions": len(net.reactions.reactions),
            "conservation": check_conservation(net),
            "dead_end_nuclides": check_dead_ends(net),
            "rate_sanity": check_rate_sanity(net),
            "xml_export": check_xml_export(case, net),
            "beta_decay_rates": check_rates(net),
        }

    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / "network_validation.json"
    path.write_text(json.dumps(report, indent=1))

    composition_report = report["composition"]
    print(f"composition: sum X = {composition_report['sum_of_mass_fractions']:.12f}")
    print(f"             X(14N) = {composition_report['x_n14']:.4e} "
          f"(proposal {composition_report['x_n14_in_proposal']:.2e})")
    print(f"             X(15N) = {composition_report['x_n15']:.4e} "
          f"(proposal {composition_report['x_n15_in_proposal']:.2e})")
    for case, record in report["networks"].items():
        conservation = record["conservation"]
        print(f"{case}: {record['species']:3d} nuclides, {record['reactions']:5d} reactions, "
              f"{len(conservation['baryon_violations'])} baryon violations, "
              f"{len(record['dead_end_nuclides'])} dead ends, "
              f"{len(record['rate_sanity']['non_finite'])} infinite rates, "
              f"{len(record['rate_sanity']['prompt_one_body_decays'])} prompt decays, "
              f"charge changes all beta: {conservation['charge_changes_are_all_beta']}")
        export = record["xml_export"]
        print(f"    XML export matches the archive: counts {export['same_counts']}, "
              f"largest rate difference {export['largest_rate_difference']}, "
              f"largest mass difference {export['largest_mass_difference']}")
        for name, values in record["beta_decay_rates"].items():
            print(f"    {name:12s} rate {values['rate_from_network']:.6e} s^-1 "
                  f"vs {values['rate_from_half_life']:.6e} from the half-life "
                  f"({values['relative_difference']:+.2%})")
    print("wrote", path)


if __name__ == "__main__":
    main()
