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
    "f17_pg_ne18": r"$^{17}$F(p,$\gamma$)",
    "ne18_bd_f18": r"$^{18}$Ne($\beta^+$)",
    "f18_pa_o15": r"$^{18}$F(p,$\alpha$)",
    "o17_pg_f18": r"$^{17}$O(p,$\gamma$)",
}

#: The CNO-II/III path, which turns out to supply most of the 15O near the
#: temperature maximum.  Drawn against 14N(p,gamma)15O so that the size of the
#: side feeding can be read straight off the figure.
SIDE_PATH = ["o16_pg_f17", "f17_pg_ne18", "ne18_bd_f18", "f18_pa_o15"]

RATIO_LABEL = r"$R_{15/14} = (^{15}$N$+^{15}$O$)/(^{14}$N$+^{14}$O$)$"

#: Every link of the closed hot CNO cycle, in cycle order.
STEADY_FLOW_CHAIN = ["c12_pg_n13", "n13_pg_o14", "o14_bd_n14",
                     "n14_pg_o15", "o15_bd_n15", "n15_pa_c12"]

#: The steady-flow plot is restricted to times at which the cycle throughput --
#: the smallest flow around the closed loop -- exceeds this fraction of its own
#: maximum.  Using the throughput rather than the reference flow alone excludes
#: the transient spike in 14N(p,gamma) that precedes the temperature maximum.
FLOW_FLOOR = 0.1


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
    """Principal cycle links above, the CNO-II/III side path below.

    The upper panel keeps the logarithmic time axis, which is the only way to
    show ten decades of evolution at once.  The lower panel is deliberately
    different: it is linear in time and restricted to the burning episode, so
    that the size of the side feeding of 15O relative to 14N(p,gamma)15O can be
    read straight off the figure.  On the logarithmic axis the whole episode is
    a single spike and the comparison is invisible.
    """
    data = read_csv(f"{run}_flows.csv")
    mask = positive(data["time"])
    shown = ["c12_pg_n13", "n14_pg_o15", "n15_pa_c12",
             "o14_bd_n14", "o15_bd_n15", "n13_pg_o14"]
    side = [r for r in SIDE_PATH if f"F_{r}" in data]

    if not side:
        # Older result files predate the CNO-II/III columns; fall back to the
        # single-panel form rather than drawing an empty axis.
        fig, axes = plt.subplots(figsize=(5.4, 3.6))
        axes = [axes]
    else:
        fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.8))
        axes = list(axes)

    for index, reaction in enumerate(shown):
        key = f"F_{reaction}"
        if key not in data:
            continue
        axes[0].plot(data["time"][mask], np.clip(data[key][mask], 1e-30, None),
                     color=PALETTE[index], label=REACTION_LABEL[reaction])
    axes[0].set_xscale("log")
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylim(1e-22, None)
    axes[0].legend(ncol=2, loc="lower left", title="closed cycle")

    if side:
        reference = data["F_n14_pg_o15"]
        # The hot phase, taken as the times at which the temperature is within a
        # factor of two of its maximum.  This is a plotting window only, chosen
        # so that the crossing of the side path over the reference flow is not
        # squeezed into a few pixels.
        inside = data["t9"] > 0.5 * data["t9"].max()
        t_lo, t_hi = data["time"][inside].min(), data["time"][inside].max()

        axes[1].plot(data["time"][inside], np.clip(reference[inside], 1e-30, None),
                     color=INK, lw=1.4,
                     label=REACTION_LABEL["n14_pg_o15"] + " (ref.)")
        for index, reaction in enumerate(side):
            axes[1].plot(data["time"][inside],
                         np.clip(data[f"F_{reaction}"][inside], 1e-30, None),
                         color=PALETTE[index], ls="--",
                         label=REACTION_LABEL[reaction])
        hottest = data["time"][int(np.argmax(data["t9"]))]
        axes[1].axvline(hottest, color=MUTED, lw=0.9, ls=(0, (2, 3)))
        axes[1].annotate(r"$T_{9,\max}$", (hottest, 1.0), xycoords=("data", "axes fraction"),
                         textcoords="offset points", xytext=(3, -10),
                         color=MUTED, fontsize=8)
        axes[1].set_xlim(t_lo, t_hi)
        axes[1].set_xlabel("time (s)")
        largest = max(reference[inside].max(),
                      max(data[f"F_{r}"][inside].max() for r in side))
        axes[1].set_ylim(1.0e-6 * largest, 5.0 * largest)
        axes[1].legend(ncol=2, loc="lower left", title="CNO-II/III side path")

    for ax in axes:
        ax.set_yscale("log")
        ax.set_ylabel(r"flow $F$ (mol g$^{-1}$ s$^{-1}$)")
    axes[0].set_title(title, loc="left", color=INK)
    fig.tight_layout()
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
    """Flow ratios for every link of the closed cycle, plus their dispersion.

    Five things this figure has to do, none of which a plot of selected
    proton-capture ratios over the whole run can do: show every link of the
    closed cycle including the two beta decays, restrict the time axis to the
    interval where the reference flow is actually large enough for the ratios
    to mean anything, mark the temperature maximum, quantify the convergence
    with a single number, and make it possible to see whether all the links are
    close to unity *at the same time*.
    """
    data = read_csv(f"{run}_flows.csv")
    time = data["time"]
    links = [f"F_{name}" for name in STEADY_FLOW_CHAIN if f"F_{name}" in data]
    if len(links) < len(STEADY_FLOW_CHAIN):
        return
    throughput = np.vstack([data[k] for k in links]).min(axis=0)
    if not np.any(throughput > 0):
        return
    # The ratios say something about circulation only while the cycle is
    # actually circulating; outside that the denominator is collapsing.
    window = throughput > FLOW_FLOOR * np.nanmax(throughput)
    if window.sum() < 3:
        return
    lo, hi = time[window][0], time[window][-1]
    span = time >= lo
    span &= time <= hi

    fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.2), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    for index, reaction in enumerate(STEADY_FLOW_CHAIN):
        key = f"Q_{reaction}"
        if key not in data:
            continue
        # The last link overlies the one before it almost exactly -- 15O decays
        # to 15N, which is consumed immediately -- so it is drawn dashed to let
        # the curve underneath show through.
        axes[0].plot(time[span], data[key][span], color=PALETTE[index],
                     linestyle=(0, (5, 2)) if index == len(STEADY_FLOW_CHAIN) - 1 else "-",
                     label=REACTION_LABEL[reaction])
    axes[0].axhline(1.0, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$Q_j = F_j / F_{^{14}\mathrm{N}(\mathrm{p},\gamma)}$")
    axes[0].set_title(title, loc="left", color=INK)
    axes[0].legend(ncol=2, loc="lower left", fontsize=7)

    dispersion = data.get("steady_flow_dispersion_dex")
    if dispersion is not None:
        axes[1].plot(time[span], dispersion[span], color=INK)
        for factor, style in ((2.0, (0, (4, 3))), (3.0, (0, (1, 2)))):
            axes[1].axhline(np.log10(factor), color=MUTED, linewidth=0.9, linestyle=style)
            # Right-aligned and offset apart, so the two labels do not collide.
            axes[1].annotate(f"factor {factor:.0f}", (0.995, np.log10(factor)),
                             xycoords=axes[1].get_yaxis_transform(),
                             textcoords="offset points", xytext=(0, 2),
                             fontsize=7, color=MUTED, va="bottom", ha="right")
        best = int(np.nanargmin(np.where(span, dispersion, np.inf)))
        axes[1].annotate(f"min $D$ = {dispersion[best]:.2f} dex",
                         (time[best], dispersion[best]), textcoords="offset points",
                         xytext=(6, 6), fontsize=8, color=INK)
    axes[1].set_ylabel(r"$D = \max_j |\log_{10} Q_j|$")
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylim(0, None)

    hottest = time[int(np.argmax(data["t9"]))]
    for ax in axes:
        ax.axvline(hottest, color=PALETTE[7], linewidth=0.9, linestyle=(0, (2, 2)))
    axes[0].annotate("$T_9$ max", (hottest, 1.0), xycoords=("data", "axes fraction"),
                     textcoords="offset points", xytext=(3, -10),
                     fontsize=8, color=PALETTE[7], va="top")
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
            marker="o", markersize=5, label="matched exponential series")
    # The two longest-exposure runs sit close together, and the trajectory point
    # and its residual arrow occupy the space directly above them, so their
    # labels are placed to the left rather than on top.
    crowded = sorted(e[0] for e in entries)[-2:] if len(entries) > 2 else []
    for x, y in entries:
        if x in crowded:
            ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                        xytext=(-7, -3), ha="right", fontsize=8, color=INK)
        else:
            ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 8),
                        ha="center", fontsize=8, color=INK)
    if trajectory:
        tx = trajectory["time_above_t9_0p2_s"]
        ty = trajectory["r_final"]
        ax.plot([tx], [ty], color=PALETTE[1], marker="D", markersize=7,
                linestyle="none", label="S1-like trajectory")
        ax.annotate(f"{ty:.2f}", (tx, ty), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8, color=PALETTE[1])
        # The residual is the gap at *equal exposure*, so it is anchored on the
        # control built to reproduce the trajectory's own exposure rather than on
        # whichever member of the series happens to come closest in value.  The
        # curve is not monotonic in exposure, so those are not the same run.
        matched = summary["runs"].get("exp_matched_exposure")
        if matched is not None:
            base = matched["r_final"]
            ax.annotate("", xy=(tx, ty), xytext=(tx, base),
                        arrowprops=dict(arrowstyle="<->", color=MUTED, linewidth=0.9))
            ax.text(tx * 1.12, (ty * base) ** 0.5, "heating phase\n+ $T$--$\\rho$ path",
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
            delta_text = f"{100.0 * delta:+.4f}\\%"
        else:
            delta_text = "--"
        # Six decimals, so that the sign of the intermediate-to-large difference
        # can be checked against the ratios themselves: at four decimals both
        # round to the same number and the quoted sign looks unsupported.
        rows.append(
            f"{name} & {limit} & {record['species']} & {record['reactions']} & "
            f"{record['r_final']:.6f} & {delta_text} \\\\"
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


def _peak_conditions(runs: dict):
    """The trajectory's temperature maximum and the density at that moment.

    Read from the matched controls, which are built from them, so that captions
    and row labels cannot drift away from the trajectory actually used.
    """
    record = runs.get("exp_matched_tau0p2") or runs.get("exp_matched_exposure")
    if record is None:
        return None, None, ""
    t9 = record["t9_0"]
    rho = record["rho_0"]
    exponent = int(np.floor(np.log10(abs(rho))))
    mantissa = rho / 10.0 ** exponent
    return t9, rho, f"{mantissa:.2f}\\times10^{{{exponent}}}"


def table_matched_control(summary: dict) -> None:
    """The matched-control series, and the decomposition it makes possible."""
    runs = summary["runs"]
    peak_t9, peak_rho, peak_rho_tex = _peak_conditions(runs)
    rows = []
    for name, record in sorted(runs.items(),
                               key=lambda kv: kv[1].get("tau", 0.0)):
        if not name.startswith("exp_matched_"):
            continue
        bold = name == "exp_matched_exposure"
        cells = [f"{record['tau']:.2f}", f"{record['time_above_t9_0p2_s']:.2f}",
                 f"{record['r_final']:.3f}", f"{record['enhancement']:.0f}"]
        if bold:
            cells = [f"\\textbf{{{c}}}" for c in cells]
        rows.append(" & ".join(cells) + " \\\\")
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
        f"maximum ($T_{{9,0}} = {peak_t9:.4f}$, $\\rho_0 = {peak_rho_tex}$ "
        "g\\,cm$^{-3}$), so that the exposure time is the only variable along "
        "the series. The row in bold is built to reproduce the trajectory's own "
        "exposure above $T_9 = 0.2$ rather than to bracket it. The trajectory "
        "result is repeated for comparison.",
        "tab:nucnetpy_matched_control",
        "$\\tau$ (s) & $t(T_9 > 0.2)$ (s) & $R^{\\mathrm{final}}_{15/14}$ & "
        "$f_{\\mathrm{enh}}$",
        rows,
        "llll",
    )

    # Labels are built from the runs themselves.  They were once written out by
    # hand, which left them quoting a previous trajectory's peak conditions after
    # the profile was replaced.
    exposure_run = "exp_matched_exposure" if "exp_matched_exposure" in runs else "exp_matched_tau7p0"
    steps = [
        (f"Peak temperature, $0.200 \\rightarrow {peak_t9:.3f}$",
         "exp_ref", "exp_peakT_only"),
        # Not "peak density": the trajectory's highest density occurs near its
        # start.  This is the density at the moment the temperature peaks.
        (f"Density at $T_{{9,\\max}}$, $1.5\\times10^{{4}} \\rightarrow {peak_rho_tex}$",
         "exp_peakT_only", "exp_matched_tau0p2"),
        (f"Exposure, ${runs['exp_matched_tau0p2']['time_above_t9_0p2_s']:.2f} \\rightarrow "
         f"{runs[exposure_run]['time_above_t9_0p2_s']:.2f}$ s",
         "exp_matched_tau0p2", exposure_run),
        ("Heating phase and $T$--$\\rho$ path", exposure_run, "traj_ref"),
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
    reversed_note = ""
    if "exp_rho_only" in runs:
        first = runs["exp_rho_only"]["r_final"] / runs["exp_ref"]["r_final"]
        second = runs["exp_matched_tau0p2"]["r_final"] / runs["exp_rho_only"]["r_final"]
        reversed_note = (
            f" Taking the first two steps in the opposite order gives "
            f"${first:.2f}$ and ${second:.2f}$ in place of the values shown."
        )
    _table(
        TABLES / "tab_decomposition.tex",
        "Sequential, order-dependent factorisation of the difference between "
        "the reference exponential model and the nova trajectory. Each row "
        "changes one property of the thermodynamic history, in the order "
        "listed; because the network is non-linear the individual factors "
        "belong to that sequence and are not a unique causal decomposition."
        + reversed_note +
        " The last row is the residual that the exponential family cannot "
        "reproduce.",
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


def table_side_flow(run: str = "traj_ref") -> None:
    """Side-flow fractions of the six cycle intermediates, over the interval in
    which all six principal links agree to within a factor of two."""
    path = RESULTS / f"{run}_flows.csv"
    if not path.exists():
        return
    columns = read_csv(path.name)
    if "fside_c12" not in columns:
        return

    dispersion = columns["steady_flow_dispersion_dex"]
    inside = np.isfinite(dispersion) & (dispersion < np.log10(2.0))
    if not inside.any():
        return
    time = columns["time"][inside]

    labels = [
        ("c12", "$^{12}\\mathrm{C}$"),
        ("n13", "$^{13}\\mathrm{N}$"),
        ("o14", "$^{14}\\mathrm{O}$"),
        ("n14", "$^{14}\\mathrm{N}$"),
        ("o15", "$^{15}\\mathrm{O}$"),
        ("n15", "$^{15}\\mathrm{N}$"),
    ]
    def scientific(value: float) -> str:
        mantissa, exponent = f"{value:.1e}".split("e")
        return f"${mantissa}\\times10^{{{int(exponent)}}}$"

    rows = []
    for key, label in labels:
        values = columns[f"fside_{key}"][inside]
        rows.append(
            f"{label} & {scientific(float(np.median(values)))} & "
            f"{scientific(float(values.max()))} \\\\"
        )

    _table(
        TABLES / "tab_side_flow.tex",
        "Side-flow fraction $f_{\\mathrm{side},i}$ of each intermediate of the "
        "closed hot-CNO cycle, over the interval "
        f"$t={time[0]:.2f}$--${time[-1]:.2f}\\ \\mathrm{{s}}$ in which all six "
        "principal links agree to within a factor of two. $f_{\\mathrm{side},i}$ "
        "is the fraction of the total flow through nuclide $i$ that is carried by "
        "reactions outside the cycle of Eq.~\\eqref{eq:closed_cycle}; values well "
        "below unity mean the nuclide is fed and drained by the cycle alone.",
        "tab:nucnetpy_side_flow",
        "Nuclide & median $f_{\\mathrm{side},i}$ & max $f_{\\mathrm{side},i}$",
        rows,
        "lll",
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
                           "Trajectory model: the six principal cycle flows")
        figure_energy("traj_ref", "fig15_traj_energy.png",
                      "Trajectory model: nuclear energy generation")
    figure_timescales("exp_ref", "fig10_exp_timescales.png",
                      "Exponential model: nuclear and thermodynamic timescales")
    figure_steady_flow("exp_ref", "fig12_exp_steady_flow.png",
                       "Exponential model: the six principal cycle flows")
    figure_energy("exp_ref", "fig14_exp_energy.png",
                  "Exponential model: nuclear energy generation")

    figure_network_size(summary, "fig16_network_size.png")
    figure_tau_series(summary, "fig17_tau_series.png")
    figure_matched_control(summary, "fig18_matched_control.png")

    table_comparison(summary)
    table_network_size(summary)
    table_tau(summary)
    table_matched_control(summary)
    table_side_flow()


if __name__ == "__main__":
    main()
