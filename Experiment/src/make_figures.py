"""Draw the figures and write the LaTeX tables for the experiment.

Reads only the CSV files and ``summary.json`` produced by ``run_experiment.py``,
so figures can be redrawn without repeating any calculation::

    python src/make_figures.py

Plotting conventions
--------------------
Each panel carries a single quantity on its y-axis; where two quantities must be
compared (temperature and density, for instance) they are drawn as stacked
panels rather than on two scales.  Series colours are assigned once, by isotope
or by reaction, and reused in every figure, so a colour always means the same
thing.  Any figure with more than one series carries a legend.  The palette is
colour-vision-safe.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
TABLES = RESULTS / "tables"

#: Colour-vision-safe categorical palette, used in a fixed order.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
           "#008300", "#4a3aa7", "#e34948"]

INK = "#1b1b1b"
MUTED = "#6b6b6b"
GRID = "#dcdcdc"

#: A colour per isotope, fixed across every figure.
ISOTOPE_COLOUR = {
    "n14": PALETTE[0], "n15": PALETTE[1], "o14": PALETTE[2], "o15": PALETTE[3],
    "c12": PALETTE[4], "c13": PALETTE[6], "o16": PALETTE[5], "f17": PALETTE[7],
}
ISOTOPE_LABEL = {
    "h1": r"$^{1}$H", "he4": r"$^{4}$He", "c12": r"$^{12}$C", "c13": r"$^{13}$C",
    "n13": r"$^{13}$N", "n14": r"$^{14}$N", "n15": r"$^{15}$N",
    "o14": r"$^{14}$O", "o15": r"$^{15}$O", "o16": r"$^{16}$O", "o17": r"$^{17}$O",
    "f17": r"$^{17}$F", "ne20": r"$^{20}$Ne",
}

REACTION_LABEL = {
    "c12_pg_n13": r"$^{12}$C(p,$\gamma$)",
    "n13_pg_o14": r"$^{13}$N(p,$\gamma$)",
    "n13_bd_c13": r"$^{13}$N($\beta^+$)",
    "c13_pg_n14": r"$^{13}$C(p,$\gamma$)",
    "n14_pg_o15": r"$^{14}$N(p,$\gamma$)",
    "o14_bd_n14": r"$^{14}$O($\beta^+$)",
    "o15_bd_n15": r"$^{15}$O($\beta^+$)",
    "n15_pa_c12": r"$^{15}$N(p,$\alpha$)",
    "n15_pg_o16": r"$^{15}$N(p,$\gamma$)",
    "o16_pg_f17": r"$^{16}$O(p,$\gamma$)",
    "o17_pa_n14": r"$^{17}$O(p,$\alpha$)",
    "f17_bd_o17": r"$^{17}$F($\beta^+$)",
}

RATIO_LABEL = r"$R_{15/14} = (^{15}$N$+^{15}$O$)/(^{14}$N$+^{14}$O$)$"


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 160,
        "savefig.dpi": 160,
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "lines.linewidth": 1.6,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "text.color": INK,
        "figure.facecolor": "white",
        "savefig.bbox": "tight",
    })


def read_csv(name: str) -> dict:
    path = RESULTS / name
    with path.open() as handle:
        header = handle.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return {key: data[:, i] for i, key in enumerate(header)}


def positive(time: np.ndarray) -> np.ndarray:
    return time > 0.0


def save(fig, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / name)
    plt.close(fig)
    print("wrote", FIGURES / name)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def figure_thermodynamics(run: str, name: str, title: str, log_time: bool) -> None:
    data = read_csv(f"{run}_evolution.csv")
    mask = positive(data["time"]) if log_time else np.ones_like(data["time"], dtype=bool)
    fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.0), sharex=True)

    axes[0].plot(data["time"][mask], data["t9"][mask], color=PALETTE[0])
    axes[0].set_ylabel(r"$T_9$")
    axes[0].set_title(title, loc="left", color=INK)

    axes[1].plot(data["time"][mask], data["rho"][mask], color=PALETTE[1])
    axes[1].set_ylabel(r"$\rho$ (g cm$^{-3}$)")
    axes[1].set_xlabel("time (s)")
    axes[1].set_yscale("log")

    if log_time:
        axes[1].set_xscale("log")
        axes[0].set_yscale("log")
    save(fig, name)


def figure_ratio(runs: list[tuple[str, str]], name: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for index, (run, label) in enumerate(runs):
        data = read_csv(f"{run}_evolution.csv")
        mask = positive(data["time"])
        ax.plot(data["time"][mask], data["R15_14"][mask],
                color=PALETTE[index], label=label)
        ax.annotate(f"{data['R15_14'][-1]:.2f}",
                    (data["time"][mask][-1], data["R15_14"][mask][-1]),
                    textcoords="offset points", xytext=(4, 0),
                    color=PALETTE[index], fontsize=8, va="center")
    initial = read_csv(f"{runs[0][0]}_evolution.csv")["R15_14"][0]
    ax.axhline(initial, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
    ax.set_xscale("log")
    ax.set_yscale("log")
    # Placed in axis coordinates along x, so the label stays inside the panel
    # whatever the time range turns out to be.
    ax.annotate(f"initial (solar) = {initial:.2e}", (0.015, initial),
                xycoords=ax.get_yaxis_transform(), textcoords="offset points",
                xytext=(0, 5), color=MUTED, fontsize=8, va="bottom", ha="left")
    ax.set_xlabel("time (s)")
    ax.set_ylabel(RATIO_LABEL)
    ax.set_title(title, loc="left", color=INK)
    if len(runs) > 1:
        ax.legend(loc="lower right")
    save(fig, name)


def figure_abundances(run: str, name: str, title: str) -> None:
    data = read_csv(f"{run}_evolution.csv")
    mask = positive(data["time"])
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for species in ["c12", "n14", "n15", "o14", "o15", "o16"]:
        key = f"X_{species}"
        if key not in data:
            continue
        values = np.clip(data[key][mask], 1e-16, None)
        ax.plot(data["time"][mask], values,
                color=ISOTOPE_COLOUR.get(species, MUTED),
                label=ISOTOPE_LABEL[species])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-12, 1e-1)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("mass fraction $X$")
    ax.set_title(title, loc="left", color=INK)
    ax.legend(ncol=3, loc="lower left")
    save(fig, name)


def figure_flows(run: str, name: str, title: str) -> None:
    data = read_csv(f"{run}_flows.csv")
    mask = positive(data["time"])
    shown = ["c12_pg_n13", "n14_pg_o15", "n15_pa_c12",
             "o14_bd_n14", "o15_bd_n15", "n13_pg_o14"]
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for index, reaction in enumerate(shown):
        key = f"F_{reaction}"
        if key not in data:
            continue
        ax.plot(data["time"][mask], np.clip(data[key][mask], 1e-30, None),
                color=PALETTE[index], label=REACTION_LABEL[reaction])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-22, None)
    ax.set_xlabel("time (s)")
    ax.set_ylabel(r"flow $F$ (mol g$^{-1}$ s$^{-1}$)")
    ax.set_title(title, loc="left", color=INK)
    ax.legend(ncol=2, loc="lower left")
    save(fig, name)


#: Timescales are plotted over this band.  A constant-temperature phase gives an
#: infinite thermodynamic timescale, and a rate that has underflowed gives an
#: infinite nuclear one; both would otherwise ruin the axis.
TAU_LIMITS = (1.0e-4, 1.0e12)


def _clipped(values: np.ndarray) -> np.ndarray:
    out = np.array(values, dtype=float)
    out[~np.isfinite(out)] = np.nan
    out[out > TAU_LIMITS[1]] = np.nan
    return out


def figure_timescales(run: str, name: str, title: str) -> None:
    data = read_csv(f"{run}_flows.csv")
    mask = positive(data["time"])
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    nuclear = ["c12_pg_n13", "n14_pg_o15", "o14_bd_n14", "o15_bd_n15"]
    for index, reaction in enumerate(nuclear):
        key = f"tau_{reaction}"
        if key not in data:
            continue
        ax.plot(data["time"][mask], _clipped(data[key][mask]),
                color=PALETTE[index], label=REACTION_LABEL[reaction])
    ax.plot(data["time"][mask], _clipped(data["tau_T9"][mask]), color=INK,
            linestyle=(0, (4, 3)), label=r"$\tau_T$")
    ax.plot(data["time"][mask], _clipped(data["tau_rho"][mask]), color=MUTED,
            linestyle=(0, (1, 2)), label=r"$\tau_\rho$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(*TAU_LIMITS)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("timescale (s)")
    ax.set_title(title, loc="left", color=INK)
    ax.legend(ncol=3, loc="lower right")
    save(fig, name)


def figure_steady_flow(run: str, name: str, title: str) -> None:
    data = read_csv(f"{run}_flows.csv")
    mask = positive(data["time"])
    chain = ["c12_pg_n13", "n13_pg_o14", "n15_pa_c12"]
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for index, reaction in enumerate(chain):
        key = f"Q_{reaction}"
        if key not in data:
            continue
        ax.plot(data["time"][mask], np.clip(data[key][mask], 1e-8, 1e8),
                color=PALETTE[index],
                label=f"{REACTION_LABEL[reaction]} / {REACTION_LABEL['n14_pg_o15']}")
    ax.axhline(1.0, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.annotate("steady flow", (0.015, 1.0),
                xycoords=ax.get_yaxis_transform(), textcoords="offset points",
                xytext=(0, 5), color=MUTED, fontsize=8, va="bottom", ha="left")
    ax.set_xlabel("time (s)")
    ax.set_ylabel(r"flow ratio $Q_{ab}$")
    ax.set_title(title, loc="left", color=INK)
    ax.legend(loc="best")
    save(fig, name)


def figure_network_size(summary: dict, name: str) -> None:
    runs = [("traj_ref", r"small, $Z\leq10$"),
            ("traj_z20", r"intermediate, $Z\leq20$"),
            ("traj_z30", r"large, $Z\leq30$")]
    runs = [(r, l) for r, l in runs if (RESULTS / f"{r}_evolution.csv").exists()]
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for index, (run, label) in enumerate(runs):
        data = read_csv(f"{run}_evolution.csv")
        mask = positive(data["time"])
        final = summary["runs"][run]["r_final"]
        ax.plot(data["time"][mask], data["R15_14"][mask], color=PALETTE[index],
                label=f"{label}  ($R_{{\\rm final}}={final:.4f}$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time (s)")
    ax.set_ylabel(RATIO_LABEL)
    ax.set_title("Trajectory model: dependence on network size", loc="left", color=INK)
    ax.legend(loc="lower right")
    save(fig, name)


def figure_tau_series(summary: dict, name: str) -> None:
    entries = []
    for run, record in summary["runs"].items():
        if run.startswith("exp_tau"):
            entries.append((record["tau"], record["r_final"]))
    if not entries:
        return
    entries.sort()
    tau = [e[0] for e in entries]
    final = [e[1] for e in entries]
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.plot(tau, final, color=PALETTE[0], marker="o", markersize=5)
    for x, y in entries:
        ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                    xytext=(0, 7), ha="center", fontsize=8, color=INK)
    ax.set_xscale("log")
    ax.set_xlabel(r"expansion timescale $\tau$ (s)")
    ax.set_ylabel(r"final $R_{15/14}$")
    ax.set_title("Exponential model: effect of the cooling rate", loc="left", color=INK)
    save(fig, name)


def figure_matched_control(summary: dict, name: str) -> None:
    """Final ratio against exposure time, at the trajectory's peak conditions.

    Everything on the curve shares one peak temperature and one peak density,
    so the only variable along it is how long the material stays hot.  The
    trajectory result is drawn as a separate point at its own exposure time;
    the gap between the two is what matching peak conditions and exposure
    cannot account for.
    """
    entries = sorted(
        (record["time_above_t9_0p2_s"], record["r_final"])
        for run, record in summary["runs"].items() if run.startswith("exp_matched_")
    )
    if not entries:
        return
    trajectory = summary["runs"].get("traj_ref")
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.plot([e[0] for e in entries], [e[1] for e in entries], color=PALETTE[0],
            marker="o", markersize=5, label="exponential, matched peak $T_9$ and $\\rho$")
    for x, y in entries:
        ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 8),
                    ha="center", fontsize=8, color=INK)
    if trajectory:
        tx = trajectory["time_above_t9_0p2_s"]
        ty = trajectory["r_final"]
        ax.plot([tx], [ty], color=PALETTE[1], marker="D", markersize=7,
                linestyle="none", label="nova trajectory")
        ax.annotate(f"{ty:.2f}", (tx, ty), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8, color=PALETTE[1])
        ax.annotate("", xy=(tx, ty), xytext=(tx, 1.9514),
                    arrowprops=dict(arrowstyle="<->", color=MUTED, linewidth=0.9))
        ax.text(tx * 1.12, (ty * 1.9514) ** 0.5, "heating phase\n+ $T$--$\\rho$ path",
                fontsize=8, color=MUTED, va="center")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"time spent above $T_9 = 0.2$ (s)")
    ax.set_ylabel(r"final $R_{15/14}$")
    ax.set_title("Exposure time at matched peak conditions", loc="left", color=INK)
    ax.legend(loc="lower right")
    save(fig, name)


def figure_energy(run: str, name: str, title: str) -> None:
    data = read_csv(f"{run}_flows.csv")
    mask = positive(data["time"])
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.plot(data["time"][mask], np.clip(data["eps_nuc_erg_g_s"][mask], 1e-6, None),
            color=PALETTE[0])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time (s)")
    ax.set_ylabel(r"$\epsilon_{\rm nuc}$ (erg g$^{-1}$ s$^{-1}$)")
    ax.set_title(title, loc="left", color=INK)
    save(fig, name)


# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------

def _table(path: Path, caption: str, label: str, header: str, rows: list[str],
           spec: str) -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    body = "\n".join(rows)
    path.write_text(
        "\\begin{table}[h!]\n\\centering\n"
        f"\\caption{{{caption}}}\n"
        f"\\begin{{tabular}}{{{spec}}}\n\\hline\n{header} \\\\\n\\hline\n"
        f"{body}\n\\hline\n\\end{{tabular}}\n"
        f"\\label{{{label}}}\n\\end{{table}}\n"
    )
    print("wrote", path)


def table_comparison(summary: dict) -> None:
    rows = []
    for run, name in [("exp_ref", "Exponential"), ("traj_ref", "Trajectory")]:
        record = summary["runs"].get(run)
        if not record:
            continue
        rows.append(
            f"{name} & {record['r_initial']:.2e} & {record['r_final']:.3f} & "
            f"{record['enhancement']:.0f} & {record['freeze_out_time_s']:.2e} & "
            f"{record['t9_peak']:.3f} \\\\"
        )
    _table(
        TABLES / "tab_comparison.tex",
        "Single-zone results for the two thermodynamic histories, computed with "
        "NucNetPy on the small ($Z\\leq10$) network.",
        "tab:nucnetpy_comparison",
        "Model & $R^{\\mathrm{initial}}_{15/14}$ & $R^{\\mathrm{final}}_{15/14}$ & "
        "$f_{\\mathrm{enh}}$ & $t_{\\mathrm{fo}}$ (s) & $T_{9,\\mathrm{max}}$",
        rows,
        "llllll",
    )


def table_network_size(summary: dict) -> None:
    reference = summary["runs"].get("traj_z30")
    rows = []
    for run, name, limit in [("traj_ref", "Small", "$Z\\leq10$"),
                             ("traj_z20", "Intermediate", "$Z\\leq20$"),
                             ("traj_z30", "Large", "$Z\\leq30$")]:
        record = summary["runs"].get(run)
        if not record:
            continue
        if reference:
            delta = (record["r_final"] - reference["r_final"]) / reference["r_final"]
            # As a percentage: the differences between network sizes are far
            # too small to read off three decimal places of a bare fraction.
            delta_text = f"{100.0 * delta:+.3f}\\%"
        else:
            delta_text = "--"
        rows.append(
            f"{name} & {limit} & {record['species']} & {record['reactions']} & "
            f"{record['r_final']:.4f} & {delta_text} \\\\"
        )
    _table(
        TABLES / "tab_network_size.tex",
        "Final CNO ratio from the trajectory model for the three network cases. "
        "$\\Delta R$ is measured against the large network.",
        "tab:nucnetpy_network_size",
        "Network & Range & Nuclides & Reactions & $R^{\\mathrm{final}}_{15/14}$ & $\\Delta R$",
        rows,
        "llllll",
    )


def table_matched_control(summary: dict) -> None:
    """The matched-control series, and the decomposition it makes possible."""
    runs = summary["runs"]
    rows = []
    for name, record in sorted(runs.items(),
                               key=lambda kv: kv[1].get("tau", 0.0)):
        if not name.startswith("exp_matched_"):
            continue
        rows.append(
            f"{record['tau']:.1f} & {record['time_above_t9_0p2_s']:.2f} & "
            f"{record['r_final']:.3f} & {record['enhancement']:.0f} \\\\"
        )
    trajectory = runs.get("traj_ref")
    if trajectory:
        rows.append("\\hline")
        rows.append(
            f"trajectory & {trajectory['time_above_t9_0p2_s']:.2f} & "
            f"{trajectory['r_final']:.3f} & {trajectory['enhancement']:.0f} \\\\"
        )
    if not rows:
        return
    _table(
        TABLES / "tab_matched_control.tex",
        "Exponential models matched to the trajectory at its temperature "
        "maximum ($T_{9,0} = 0.4481$, $\\rho_0 = 4.07\\times10^{3}$ "
        "g\\,cm$^{-3}$), so that the exposure time is the only variable along "
        "the series. The trajectory result is repeated for comparison.",
        "tab:nucnetpy_matched_control",
        "$\\tau$ (s) & $t(T_9 > 0.2)$ (s) & $R^{\\mathrm{final}}_{15/14}$ & "
        "$f_{\\mathrm{enh}}$",
        rows,
        "llll",
    )

    steps = [
        ("Peak temperature, $0.200 \\rightarrow 0.448$", "exp_ref", "exp_peakT_only"),
        ("Peak density, $1.5\\times10^{4} \\rightarrow 4.07\\times10^{3}$",
         "exp_peakT_only", "exp_matched_tau0p2"),
        ("Exposure, $0.48 \\rightarrow 16.9$ s", "exp_matched_tau0p2", "exp_matched_tau7p0"),
        ("Heating phase and $T$--$\\rho$ path", "exp_matched_tau7p0", "traj_ref"),
    ]
    rows, product = [], 1.0
    for label, before, after in steps:
        if before not in runs or after not in runs:
            return
        factor = runs[after]["r_final"] / runs[before]["r_final"]
        product *= factor
        rows.append(
            f"{label} & {runs[before]['r_final']:.3f} & "
            f"{runs[after]['r_final']:.3f} & {factor:.2f} \\\\"
        )
    rows.append("\\hline")
    rows.append(
        f"Combined & {runs['exp_ref']['r_final']:.3f} & "
        f"{runs['traj_ref']['r_final']:.3f} & {product:.2f} \\\\"
    )
    _table(
        TABLES / "tab_decomposition.tex",
        "Decomposition of the difference between the reference exponential "
        "model and the nova trajectory into single-variable steps. Each row "
        "changes one property of the thermodynamic history; the last row is "
        "the residual that the exponential family cannot reproduce.",
        "tab:nucnetpy_decomposition",
        "Change & $R$ before & $R$ after & Factor",
        rows,
        "llll",
    )


def table_tau(summary: dict) -> None:
    entries = sorted(
        (record["tau"], record) for run, record in summary["runs"].items()
        if run.startswith("exp_tau")
    )
    rows = [
        f"{tau:.2f} & {record['r_final']:.3f} & {record['enhancement']:.0f} & "
        f"{record['freeze_out_time_s']:.2e} \\\\"
        for tau, record in entries
    ]
    if not rows:
        return
    _table(
        TABLES / "tab_tau.tex",
        "Exponential model: effect of the expansion timescale $\\tau$ on the final "
        "CNO ratio and the freeze-out time.",
        "tab:nucnetpy_tau",
        "$\\tau$ (s) & $R^{\\mathrm{final}}_{15/14}$ & $f_{\\mathrm{enh}}$ & "
        "$t_{\\mathrm{fo}}$ (s)",
        rows,
        "llll",
    )


def main() -> None:
    style()
    summary = json.loads((RESULTS / "summary.json").read_text())

    figure_thermodynamics("exp_ref", "fig01_exp_thermodynamics.png",
                          "Exponential model: thermodynamic history", log_time=True)
    figure_ratio([("exp_ref", "exponential")], "fig02_exp_ratio.png",
                 "Exponential model: CNO diagnostic ratio")
    figure_abundances("exp_ref", "fig03_exp_abundances.png",
                      "Exponential model: CNO isotopes")

    if (RESULTS / "traj_ref_evolution.csv").exists():
        figure_thermodynamics("traj_ref", "fig04_traj_thermodynamics.png",
                              "Trajectory model: thermodynamic history", log_time=True)
        figure_ratio([("traj_ref", "trajectory")], "fig05_traj_ratio.png",
                     "Trajectory model: CNO diagnostic ratio")
        figure_abundances("traj_ref", "fig06_traj_abundances.png",
                          "Trajectory model: CNO isotopes")
        figure_ratio([("exp_ref", "exponential"), ("traj_ref", "trajectory")],
                     "fig07_ratio_comparison.png",
                     "Effect of the thermodynamic history on $R_{15/14}$")

    figure_flows("exp_ref", "fig08_exp_flows.png", "Exponential model: hot-CNO flows")
    if (RESULTS / "traj_ref_flows.csv").exists():
        figure_flows("traj_ref", "fig09_traj_flows.png", "Trajectory model: hot-CNO flows")
        figure_timescales("traj_ref", "fig11_traj_timescales.png",
                          "Trajectory model: nuclear and thermodynamic timescales")
        figure_steady_flow("traj_ref", "fig13_traj_steady_flow.png",
                           "Trajectory model: test for steady flow")
        figure_energy("traj_ref", "fig15_traj_energy.png",
                      "Trajectory model: nuclear energy generation")
    figure_timescales("exp_ref", "fig10_exp_timescales.png",
                      "Exponential model: nuclear and thermodynamic timescales")
    figure_steady_flow("exp_ref", "fig12_exp_steady_flow.png",
                       "Exponential model: test for steady flow")
    figure_energy("exp_ref", "fig14_exp_energy.png",
                  "Exponential model: nuclear energy generation")

    figure_network_size(summary, "fig16_network_size.png")
    figure_tau_series(summary, "fig17_tau_series.png")
    figure_matched_control(summary, "fig18_matched_control.png")

    table_comparison(summary)
    table_network_size(summary)
    table_tau(summary)
    table_matched_control(summary)


if __name__ == "__main__":
    main()
