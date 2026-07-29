"""Run the nova hydrogen-burning experiment described in the proposal.

The driver carries out the computational plan of the methodology chapter:

1.  the reference exponential single-zone calculation,
2.  the reference trajectory-based single-zone calculation,
3.  exponential calculations for a range of expansion timescales ``tau``,
4.  the trajectory calculation repeated for three network sizes,
5.  the CNO reaction flows, nuclear and thermodynamic timescales, freeze-out
    times and steady-flow diagnostics for the two reference runs.

All the physics is done by NucNetPy: the network comes from
:class:`nucnetpy.Network`, the integration from :func:`nucnetpy.evolve_zone`
with its analytic Jacobian, the flows from :meth:`nucnetpy.Reaction.flux` and
the energy generation from
:func:`nucnetpy.analysis.nuclear_energy_generation_rate`.

Usage::

    python src/run_experiment.py              # every run
    python src/run_experiment.py exp_ref      # one run by name
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np

from nucnetpy import Zone, evolve_zone

import composition
import diagnostics
import thermodynamics as th
from network_io import load_network

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

#: Solver settings.  ``rtol`` was tightened until the final diagnostic ratio
#: stopped moving in the sixth significant figure.  ``atol`` is deliberately not
#: pushed further down: the abundances that matter here are never smaller than
#: about 1e-12, and asking BDF to resolve components at 1e-30 collapses the step
#: size across the thermonuclear runaway without changing any answer.
RTOL = 1.0e-8
ATOL = 1.0e-25


@dataclass
class Run:
    """One single-zone calculation."""

    name: str
    network: str
    description: str
    thermo: Callable
    t_start: float
    t_end: float
    n_times: int = 400
    detailed: bool = False
    extra: Dict[str, float] = field(default_factory=dict)
    #: Optional ``(t_first, t_last)`` window given extra, linearly spaced output
    #: points, so that a short burning episode is resolved in the output as well
    #: as in the integration.
    dense_window: Optional[tuple] = None

    def times(self) -> np.ndarray:
        grid = np.concatenate(
            ([0.0], np.logspace(np.log10(self.t_start), np.log10(self.t_end), self.n_times))
        )
        if self.dense_window is not None:
            first, last = self.dense_window
            grid = np.concatenate((grid, np.linspace(first, last, 200)))
        return np.unique(grid)


def build_runs() -> List[Run]:
    trajectory = th.nova_trajectory()
    # The trajectory file sets its own span.  Integrating past its last row
    # would mean holding the final temperature and density fixed for ever,
    # which is an extrapolation rather than a result; the composition has
    # frozen long before then in any case.
    traj_end = float(trajectory.time[-1])
    traj_peak = th.peak_time(trajectory)
    # Extra output points across the burning episode, which lasts of order a
    # minute and would otherwise fall between two points of a logarithmic grid
    # spanning five decades.
    traj_window = (0.5 * traj_peak, 5.0 * traj_peak)
    runs: List[Run] = [
        Run(
            name="exp_ref",
            network="nova_z10",
            description="Reference exponential model, T9_0=0.20, rho_0=1.5e4, tau=0.2 s",
            thermo=th.exponential_thermo(t9_0=0.20, rho_0=1.5e4, tau=0.2),
            t_start=1.0e-9,
            t_end=100.0,
            detailed=True,
            extra={"t9_0": 0.20, "rho_0": 1.5e4, "tau": 0.2},
        ),
        Run(
            name="traj_ref",
            network="nova_z10",
            description="Reference trajectory model, small network (Z <= 10)",
            thermo=trajectory.thermo,
            t_start=1.0e-5,
            t_end=traj_end,
            detailed=True,
            dense_window=traj_window,
        ),
        Run(
            name="traj_z20",
            network="nova_z20",
            description="Trajectory model, intermediate network (Z <= 20)",
            thermo=trajectory.thermo,
            t_start=1.0e-5,
            t_end=traj_end,
            dense_window=traj_window,
        ),
        Run(
            name="traj_z30",
            network="nova_z30",
            description="Trajectory model, large network (Z <= 30)",
            thermo=trajectory.thermo,
            t_start=1.0e-5,
            t_end=traj_end,
            dense_window=traj_window,
        ),
    ]

    for tau in (0.05, 0.10, 0.20, 0.50, 1.00):
        runs.append(
            Run(
                name=f"exp_tau{tau:.2f}".replace(".", "p"),
                network="nova_z10",
                description=f"Exponential model with tau = {tau:.2f} s",
                thermo=th.exponential_thermo(t9_0=0.20, rho_0=1.5e4, tau=tau),
                t_start=1.0e-9,
                t_end=100.0,
                extra={"t9_0": 0.20, "rho_0": 1.5e4, "tau": tau},
            )
        )
    return runs


def peak_temperature(thermo: Callable, t_end: float) -> float:
    """Highest ``T9`` the history reaches, on a grid fine enough to catch a spike."""
    grid = np.concatenate(([0.0], np.logspace(-9.0, np.log10(t_end), 200000)))
    return float(np.max([thermo(t)[0] for t in grid]))


def execute(run: Run) -> dict:
    """Integrate one run and write its output files.  Returns its summary."""
    net = load_network(run.network)
    y0 = composition.solar_abundances(species_names=set(net.species))
    zone = Zone(label=(run.name, "0", "0"), abundances=dict(y0))
    net.add_zone(zone)

    times = run.times()
    started = time.time()
    result = evolve_zone(net, zone, times, thermo=run.thermo, method="bdf",
                         rtol=RTOL, atol=ATOL)
    elapsed = time.time() - started

    index = {name: i for i, name in enumerate(result.species)}
    history = {name: result.y[:, i] for name, i in index.items()}
    t9 = np.array([run.thermo(t)[0] for t in result.time])
    rho = np.array([run.thermo(t)[1] for t in result.time])

    ratio = diagnostics.ratio_r15_14(history)
    summary = {
        "name": run.name,
        "network": run.network,
        "description": run.description,
        "species": len(net.species),
        "reactions": len(net.reactions.reactions),
        "solver_success": bool(result.success),
        "solver_message": result.message,
        "wall_time_s": round(elapsed, 1),
        "t_end": float(result.time[-1]),
        "reached_t_end": bool(result.time[-1] >= times[-1] * (1.0 - 1.0e-9)),
        "t9_initial": float(t9[0]),
        # The peak seen by the output grid can undershoot a short spike, so the
        # history itself is resampled finely to report the temperature the
        # integration actually experienced.  The grid has to be logarithmic:
        # a linear one over 3e7 s steps straight over a 100 s temperature spike.
        "t9_peak": peak_temperature(run.thermo, times[-1]),
        "t9_peak_on_output_grid": float(np.max(t9)),
        "rho_initial": float(rho[0]),
        "r_initial": float(ratio[0]),
        "r_final": float(ratio[-1]),
        "r_max": float(np.nanmax(ratio)),
        "enhancement": float(ratio[-1] / ratio[0]),
        "freeze_out_time_s": diagnostics.freeze_out_time(result.time, ratio),
        "strongest_change_time_s": diagnostics.time_of_strongest_change(result.time, ratio),
        "xsum_initial": float(sum(net.species[s].a * history[s][0] for s in index)),
        "xsum_final": float(sum(net.species[s].a * history[s][-1] for s in index)),
        "composition_truncation": composition.truncation_report(set(net.species)),
        **run.extra,
    }

    tracked = [s for s in diagnostics.TRACKED_SPECIES if s in index]
    columns = {"time": result.time, "t9": t9, "rho": rho, "R15_14": ratio}
    for name in tracked:
        columns[f"Y_{name}"] = history[name]
        columns[f"X_{name}"] = history[name] * net.species[name].a
    _write_csv(RESULTS / f"{run.name}_evolution.csv", columns)

    summary["final_mass_fractions"] = {
        name: float(history[name][-1] * net.species[name].a)
        for name in sorted(index, key=lambda s: (net.species[s].z, net.species[s].a))
        if history[name][-1] * net.species[name].a > 1.0e-12
    }

    if run.detailed:
        flow = diagnostics.flow_history(net, result.time, t9, rho, history, list(index))
        flow_columns = {"time": result.time, "t9": t9, "rho": rho,
                        "eps_nuc_erg_g_s": flow.energy_generation}
        for name, values in flow.flows.items():
            flow_columns[f"F_{name}"] = values
        for name, values in flow.nuclear_timescales.items():
            flow_columns[f"tau_{name}"] = values
        tau_t9, tau_rho = _thermodynamic_timescales(run.thermo, result.time)
        flow_columns["tau_T9"] = tau_t9
        flow_columns["tau_rho"] = tau_rho
        for name, values in diagnostics.steady_flow_ratios(flow).items():
            flow_columns[f"Q_{name}"] = values
        _write_csv(RESULTS / f"{run.name}_flows.csv", flow_columns)

        summary["peak_energy_generation_erg_g_s"] = float(np.nanmax(flow.energy_generation))
        summary["dominant_flow_at_peak"] = _dominant_flow(flow, t9)

    return summary


def _thermodynamic_timescales(thermo: Callable, times: np.ndarray):
    """``|T9 / dT9dt|`` and ``|rho / drhodt|`` by central differences."""
    t9 = np.array([thermo(t)[0] for t in times])
    rho = np.array([thermo(t)[1] for t in times])
    with np.errstate(divide="ignore", invalid="ignore"):
        tau_t9 = np.abs(t9 / np.gradient(t9, times))
        tau_rho = np.abs(rho / np.gradient(rho, times))
    return np.nan_to_num(tau_t9, nan=np.inf), np.nan_to_num(tau_rho, nan=np.inf)


def _dominant_flow(flow: diagnostics.FlowHistory, t9: np.ndarray) -> str:
    """Name of the strongest CNO flow at the hottest point of the run."""
    hottest = int(np.argmax(t9))
    ranked = sorted(flow.flows.items(), key=lambda item: item[1][hottest], reverse=True)
    return ranked[0][0] if ranked else ""


def _write_csv(path: Path, columns: Dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    data = np.column_stack([np.asarray(columns[n], dtype=float) for n in names])
    with path.open("w") as handle:
        handle.write(",".join(names) + "\n")
        for row in data:
            handle.write(",".join(f"{v:.8e}" for v in row) + "\n")


def main(selected: Optional[List[str]] = None) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    runs = build_runs()
    if selected:
        runs = [r for r in runs if r.name in selected]

    for run in runs:
        print(f"--- {run.name}: {run.description}", flush=True)
        summary = execute(run)
        # One file per run, so that runs can be executed as separate processes.
        (RESULTS / f"summary_{run.name}.json").write_text(json.dumps(summary, indent=1))
        print(
            f"    R_initial={summary['r_initial']:.4e}  "
            f"R_final={summary['r_final']:.4e}  "
            f"f_enh={summary['enhancement']:.1f}  "
            f"({summary['wall_time_s']:.0f} s, success={summary['solver_success']})",
            flush=True,
        )
    merge_summaries()


def merge_summaries() -> Path:
    """Collect the per-run summaries into ``results/summary.json``.

    The peak temperature is recomputed here from each run's own thermodynamic
    history, so that merging always yields a consistent summary even if a run
    was carried out before that quantity was defined as it is now.
    """
    definitions = {run.name: run for run in build_runs()}
    runs = {}
    for path in sorted(RESULTS.glob("summary_*.json")):
        summary = json.loads(path.read_text())
        run = definitions.get(summary["name"])
        if run is not None:
            summary["t9_peak"] = peak_temperature(run.thermo, summary["t_end"])
        runs[summary["name"]] = summary
    summary_path = RESULTS / "summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "solver": {"method": "bdf", "rtol": RTOL, "atol": ATOL,
                           "jacobian": "analytic"},
                "composition": "Bergemann, Lodders & Palme (2025) solar",
                "rates": "JINA ReacLib snapshot 20180319default2",
                "runs": runs,
            },
            indent=1,
        )
    )
    print("wrote", summary_path)
    return summary_path


if __name__ == "__main__":
    main(sys.argv[1:] or None)
