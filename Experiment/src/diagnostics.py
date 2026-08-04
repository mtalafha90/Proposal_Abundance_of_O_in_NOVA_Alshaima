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
    # The CNO-II/III path.  It turns out to dominate the supply of 15O near the
    # temperature maximum, so its three links are followed alongside the cycle
    # proper: 16O(p,g)17F is already in the list above.
    ("f17_pg_ne18", ["f17", "h1"], ["ne18"]),
    ("ne18_bd_f18", ["ne18"], ["f18"]),
    ("f18_pa_o15", ["f18", "h1"], ["o15", "he4"]),
    ("o17_pg_f18", ["o17", "h1"], ["f18"]),
]

#: Every link of the closed hot CNO cycle,
#:
#:   12C(p,g)13N(p,g)14O(b+)14N(p,g)15O(b+)15N(p,a)12C
#:
#: in cycle order.  A steady-flow test has to include all of them: in the
#: beta-limited regime the two decays are the slow links, so a test built only
#: from the proton captures cannot show whether the cycle is circulating at a
#: uniform rate.
STEADY_FLOW_CHAIN = [
    "c12_pg_n13", "n13_pg_o14", "o14_bd_n14",
    "n14_pg_o15", "o15_bd_n15", "n15_pa_c12",
]

#: Isotopes whose evolution is written out for every run.
TRACKED_SPECIES = [
    "h1", "he4", "c12", "c13", "n13", "n14", "n15",
    "o14", "o15", "o16", "o17", "o18", "f17", "f18", "ne18", "ne20",
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
    """Flow ratios ``Q_j = F_j / F_ref`` for every link of the closed cycle.

    The reference is :sup:`14`\\ N(p,gamma):sup:`15`\\ O, which is the limiting
    step of the cold cycle; the ratios test whether it stays so.  Values near
    one for *all* links at the *same* time mean the cycle is circulating at a
    uniform rate.
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


def steady_flow_dispersion(history: FlowHistory) -> np.ndarray:
    """``D = max_j |log10 Q_j|`` over the closed cycle: one number per time.

    ``D = 0`` is exact steady flow; ``D = log10(2)`` means every link is within
    a factor of two of the reference.  Taking the maximum over links, rather
    than an average, is what makes this a test of the whole cycle: it cannot be
    made small by a subset of the flows agreeing among themselves.
    """
    ratios = steady_flow_ratios(history)
    if not ratios:
        return np.full(len(history.time), np.nan)
    stack = np.vstack([ratios[name] for name in STEADY_FLOW_CHAIN if name in ratios])
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.nanmax(np.abs(np.log10(stack)), axis=0)


def cycle_throughput(history: FlowHistory) -> np.ndarray:
    """The rate at which the closed cycle actually circulates.

    A cycle turns over no faster than its slowest link, so the minimum flow
    around the loop is the throughput.  This is the quantity that says whether
    the cycle is running, and it is a better window criterion than the
    reference flow alone: ``14N(p,gamma)`` spikes before the temperature
    maximum while its target is being consumed, which is a transient rather
    than circulation.
    """
    stack = np.vstack([history.flows[name] for name in STEADY_FLOW_CHAIN
                       if name in history.flows])
    return np.min(stack, axis=0)


def steady_flow_report(history: FlowHistory, t9: np.ndarray,
                       floor: float = 0.1) -> Dict[str, object]:
    """Summarise how closely, and for how long, the cycle approaches steady flow.

    The analysis is restricted to times at which the cycle throughput exceeds
    ``floor`` times its own maximum.  Outside that window the flow ratios are
    dominated by a collapsing denominator and say nothing about circulation.
    """
    if "n14_pg_o15" not in history.flows:
        return {}
    throughput = cycle_throughput(history)
    if not np.any(throughput > 0):
        return {}
    dispersion = steady_flow_dispersion(history)
    window = throughput > floor * np.nanmax(throughput)
    time = history.time
    best = int(np.nanargmin(np.where(window, dispersion, np.inf)))
    out = {
        "window_criterion": f"cycle throughput > {floor:g} of its maximum",
        "window_start_s": float(time[window][0]),
        "window_end_s": float(time[window][-1]),
        "min_dispersion_dex": float(dispersion[best]),
        "min_dispersion_as_factor": float(10.0 ** dispersion[best]),
        "time_of_min_dispersion_s": float(time[best]),
        "t9_at_min_dispersion": float(t9[best]),
    }
    for factor in (1.2, 1.5, 2.0, 3.0):
        inside = window & (dispersion < np.log10(factor))
        key = f"interval_within_factor_{factor:g}_s".replace(".", "p")
        out[key] = float(time[inside][-1] - time[inside][0]) if inside.any() else 0.0
    return out


#: The intermediate nuclides of the closed cycle, each with the two principal
#: links that produce and destroy it.
CYCLE_NUCLIDES = {
    "c12": ("n15_pa_c12", "c12_pg_n13"),
    "n13": ("c12_pg_n13", "n13_pg_o14"),
    "o14": ("n13_pg_o14", "o14_bd_n14"),
    "n14": ("o14_bd_n14", "n14_pg_o15"),
    "o15": ("n14_pg_o15", "o15_bd_n15"),
    "n15": ("o15_bd_n15", "n15_pa_c12"),
}


def side_flow_fractions(net: Network, time: np.ndarray, t9: np.ndarray,
                        rho: np.ndarray, abundance_history: Dict[str, np.ndarray],
                        species_order: Sequence[str]) -> Dict[str, np.ndarray]:
    """Fraction of each cycle nuclide's turnover carried by non-cycle channels.

    For nuclide ``i``,

        f_side,i = sum |F_side,i| / (sum |F_principal,i| + sum |F_side,i|),

    where the principal flows are the two links of the closed cycle that
    produce and destroy it, and the side flows are every other reaction in the
    network that touches it, weighted by its stoichiometric coefficient.

    Equal flows around the six principal links do not by themselves establish
    steady flow of the cycle: if a nuclide is also being fed or drained from
    outside, the loop is not closed.  This quantity is what decides whether the
    stronger statement is warranted, and it has to be computed over the whole
    network rather than the tracked subset.
    """
    principal = {}
    for name, reactants, products in CNO_REACTIONS:
        reaction = find_reaction(net, reactants, products)
        if reaction is not None:
            principal[name] = reaction

    # Every reaction that touches each cycle nuclide, with its stoichiometry.
    touching = {nuc: [] for nuc in CYCLE_NUCLIDES}
    for reaction in net.reactions.reactions:
        stoichiometry = reaction.stoichiometry()
        for nuc in CYCLE_NUCLIDES:
            nu = stoichiometry.get(nuc, 0)
            if nu:
                touching[nuc].append((reaction, abs(nu)))

    out = {nuc: np.zeros(len(time)) for nuc in CYCLE_NUCLIDES}
    for step in range(len(time)):
        y = {name: abundance_history[name][step] for name in species_order}
        for nuc, (into, outof) in CYCLE_NUCLIDES.items():
            principal_ids = {id(principal[k]) for k in (into, outof) if k in principal}
            main = side = 0.0
            for reaction, nu in touching[nuc]:
                flux = nu * reaction.flux(y, t9=t9[step], rho=rho[step])
                if id(reaction) in principal_ids:
                    main += flux
                else:
                    side += flux
            total = main + side
            out[nuc][step] = side / total if total > 0.0 else 0.0
    return out


def flow_history_from_columns(columns: Dict[str, np.ndarray]) -> FlowHistory:
    """Rebuild a :class:`FlowHistory` from a stored ``*_flows.csv``.

    The steady-flow diagnostics are pure post-processing of the flows, so they
    can be recomputed from stored output without repeating an integration.
    """
    flows = {key[2:]: values for key, values in columns.items() if key.startswith("F_")}
    timescales = {key[4:]: values for key, values in columns.items() if key.startswith("tau_")}
    return FlowHistory(columns["time"], flows, timescales,
                       columns.get("eps_nuc_erg_g_s", np.zeros(len(columns["time"]))))
