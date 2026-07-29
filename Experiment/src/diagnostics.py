"""Diagnostic quantities for the nova hydrogen-burning calculations.

Everything here follows the definitions in the methodology chapter of the
proposal: the CNO ratio, its freeze-out time, the reaction flows around the
hot-CNO cycle, the nuclear and thermodynamic timescales, and the flow ratios
that test for quasi-steady flow.  The flows and rates themselves are evaluated
by NucNetPy; this module only assembles them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence

import numpy as np

from nucnetpy import Network
from nucnetpy.analysis import nuclear_energy_generation_rate

from network_io import find_reaction

#: The reactions of the hot-CNO cycle that the proposal asks to be followed.
#: Each entry is (short name, reactants, products).  The short names avoid
#: commas and brackets so that they can be used as CSV column names; the
#: readable forms live in ``make_figures.py``.
CNO_REACTIONS: List[tuple] = [
    ("c12_pg_n13", ["c12", "h1"], ["n13"]),
    ("n13_pg_o14", ["n13", "h1"], ["o14"]),
    ("n13_bd_c13", ["n13"], ["c13"]),
    ("c13_pg_n14", ["c13", "h1"], ["n14"]),
    ("n14_pg_o15", ["n14", "h1"], ["o15"]),
    ("o14_bd_n14", ["o14"], ["n14"]),
    ("o15_bd_n15", ["o15"], ["n15"]),
    ("n15_pa_c12", ["n15", "h1"], ["c12", "he4"]),
    ("n15_pg_o16", ["n15", "h1"], ["o16"]),
    ("o16_pg_f17", ["o16", "h1"], ["f17"]),
    ("o17_pa_n14", ["o17", "h1"], ["n14", "he4"]),
    ("f17_bd_o17", ["f17"], ["o17"]),
]

#: The four flows whose equality defines steady flow around the cold CNO cycle.
STEADY_FLOW_CHAIN = ["c12_pg_n13", "n13_pg_o14", "n14_pg_o15", "n15_pa_c12"]

#: Isotopes whose evolution is written out for every run.
TRACKED_SPECIES = [
    "h1", "he4", "c12", "c13", "n13", "n14", "n15",
    "o14", "o15", "o16", "o17", "o18", "f17", "f18", "ne20",
]


def ratio_r15_14(abundances: Dict[str, np.ndarray]) -> np.ndarray:
    """The diagnostic ratio ``(Y15N + Y15O) / (Y14N + Y14O)``."""
    numerator = abundances["n15"] + abundances["o15"]
    denominator = abundances["n14"] + abundances["o14"]
    return numerator / np.where(denominator > 0.0, denominator, np.nan)


def freeze_out_time(time: np.ndarray, ratio: np.ndarray, tolerance: float = 0.01) -> float:
    """Earliest time after which ``R`` stays within ``tolerance`` of its final value.

    This is the definition given in the methodology chapter: the first time
    ``t_fo`` such that ``|R(t_final) - R(t)| / R(t_final) < tolerance`` holds for
    every later time.
    """
    final = ratio[-1]
    if not np.isfinite(final) or final == 0.0:
        return float("nan")
    within = np.abs((final - ratio) / final) < tolerance
    # Walk back from the end to find where the band is entered for good.
    index = len(within) - 1
    while index > 0 and within[index - 1]:
        index -= 1
    return float(time[index])


def time_of_strongest_change(time: np.ndarray, ratio: np.ndarray) -> float:
    """Time at which ``d ln R / d ln t`` is largest."""
    positive = time > 0.0
    log_t = np.log(time[positive])
    log_r = np.log(np.clip(ratio[positive], 1.0e-300, None))
    if len(log_t) < 3:
        return float("nan")
    slope = np.gradient(log_r, log_t)
    return float(time[positive][int(np.nanargmax(slope))])


@dataclass
class FlowHistory:
    """Reaction flows and timescales sampled along a trajectory."""

    time: np.ndarray
    flows: Dict[str, np.ndarray]
    nuclear_timescales: Dict[str, np.ndarray]
    energy_generation: np.ndarray


def flow_history(
    net: Network,
    time: np.ndarray,
    t9: np.ndarray,
    rho: np.ndarray,
    abundance_history: Dict[str, np.ndarray],
    species_order: Sequence[str],
) -> FlowHistory:
    """Evaluate the hot-CNO flows and timescales at every output time.

    The flow of each reaction is NucNetPy's own :meth:`Reaction.flux`, which is
    ``rho N_A<sigma v> Y_i Y_p`` for a proton capture and ``lambda Y_i`` for a
    beta decay, exactly the two expressions of the methodology chapter.
    """
    reactions = {}
    for name, reactants, products in CNO_REACTIONS:
        reaction = find_reaction(net, reactants, products)
        if reaction is not None:
            reactions[name] = reaction

    flows = {name: np.zeros(len(time)) for name in reactions}
    timescales = {name: np.full(len(time), np.inf) for name in reactions}
    energy = np.zeros(len(time))

    for step in range(len(time)):
        y = {name: abundance_history[name][step] for name in species_order}
        y_p = y.get("h1", 0.0)
        for name, reaction in reactions.items():
            rate = reaction.rate(t9[step])
            flows[name][step] = reaction.flux(y, t9=t9[step], rho=rho[step])
            if reaction.reactant_order == 2:
                # Proton-capture timescale 1 / (rho Y_p N_A<sigma v>).
                denominator = rho[step] * y_p * rate
            else:
                # Beta-decay timescale 1 / lambda.
                denominator = rate
            timescales[name][step] = 1.0 / denominator if denominator > 0.0 else np.inf
        energy[step] = nuclear_energy_generation_rate(
            net, 0, t9=t9[step], rho=rho[step], abundances=y
        )

    return FlowHistory(time, flows, timescales, energy)


def steady_flow_ratios(history: FlowHistory) -> Dict[str, np.ndarray]:
    """Flow ratios ``Q_ab = F_a / F_b`` along the cold-CNO chain.

    Each flow is divided by the flow of :sup:`14`\\ N(p,gamma):sup:`15`\\ O, the
    classical bottleneck of the cycle.  Values near one mean the cycle is close
    to steady flow.
    """
    reference = history.flows.get("n14_pg_o15")
    if reference is None:
        return {}
    out = {}
    for name in STEADY_FLOW_CHAIN:
        if name in history.flows:
            with np.errstate(divide="ignore", invalid="ignore"):
                out[name] = np.where(reference > 0.0, history.flows[name] / reference, np.nan)
    return out
