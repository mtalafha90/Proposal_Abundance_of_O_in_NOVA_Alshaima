"""The two thermodynamic histories used in the study.

Both are represented as a NucNetPy :class:`nucnetpy.hydro.Trajectory`, whose
``thermo`` method is exactly the ``f(t, y) -> (T9, rho)`` callable that
:func:`nucnetpy.evolve_zone` expects.

Exponential model
-----------------
The analytic prescription of the proposal,

    rho(t) = rho_0 exp(-t / tau)
    T9(t)  = T9_0  exp(-t / 3 tau)

Note that NucNetPy's own :func:`nucnetpy.hydro.exponential_expansion` cools the
temperature on the same timescale as the density, so the trajectory is built
here from :meth:`Trajectory.from_columns` instead, to follow the proposal.

Trajectory model
----------------
The nova trajectory is a tabulated temperature-density history read with
:func:`nucnetpy.read_trajectory`.  If no measured trajectory file is supplied,
:func:`write_reference_trajectory` generates the reference nova history used in
this study, which reproduces the behaviour described in the proposal: a long
quiescent phase near ``T9 = 0.091`` and ``rho = 2.21e4 g/cm^3``, a
thermonuclear runaway that peaks at ``T9 = 0.447`` about 100 s in, and
adiabatic expansion and cooling afterwards, followed out to one year.

The runaway is modelled as a steep power-law rise in temperature at very nearly
constant density (the envelope is degenerate, so it cannot expand while it
heats).  After the peak the envelope expands homologously, its radius growing
linearly with time on a timescale ``tau_exp``, so that

    rho = rho_peak (1 + (t - t_peak) / tau_exp)^-3
    T9  = T9_peak  (rho / rho_peak)^(2/3)

the second relation being adiabatic expansion of an ideal monatomic gas.  This
is a reconstruction, not a hydrodynamic result: the file is deliberately kept
as a separate, swappable data file so that a trajectory extracted from a real
nova model can be dropped in without touching any of the analysis.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from nucnetpy import Trajectory, read_trajectory

ROOT = Path(__file__).resolve().parent.parent
TRAJECTORY_DIR = ROOT / "data" / "trajectories"
REFERENCE_TRAJECTORY = TRAJECTORY_DIR / "nova_reference.txt"

#: Temperature floor applied to both histories.
#:
#: ReacLib's seven-parameter fits are made for ``T9 >= 0.01`` and must not be
#: extrapolated below it: the ``a1/T9`` and ``a2/T9^(1/3)`` terms turn some
#: neutron-induced rates into 10^24 s^-1 or into infinity by ``T9 = 0.001``.
#: The floor costs nothing physically.  Charged-particle reactions are already
#: dead at ``T9 = 0.01`` -- the Coulomb barrier sees to that -- and beta decays,
#: which are what still matters during the long cooling tail, do not depend on
#: temperature at all.
T9_FLOOR = 1.0e-2

#: Density floor, low enough to be irrelevant to every two-body flow and high
#: enough to keep logarithms finite.
RHO_FLOOR = 1.0e-30

# Parameters of the reference nova trajectory.
NOVA_T9_INITIAL = 0.091
NOVA_RHO_INITIAL = 2.21e4
NOVA_T9_PEAK = 0.447
NOVA_T_PEAK = 100.0
NOVA_RUNAWAY_INDEX = 12.0
NOVA_EXPANSION_TIME = 30.0
NOVA_T_END = 3.15e7


def exponential_trajectory(
    t9_0: float = 0.20,
    rho_0: float = 1.5e4,
    tau: float = 0.2,
    t_end: float = 100.0,
    n: int = 4000,
) -> Trajectory:
    """Return the exponential expansion history of the proposal."""
    time = np.concatenate(([0.0], np.logspace(-9.0, np.log10(t_end), n - 1)))
    t9 = np.maximum(t9_0 * np.exp(-time / (3.0 * tau)), T9_FLOOR)
    rho = np.maximum(rho_0 * np.exp(-time / tau), RHO_FLOOR)
    return Trajectory.from_columns(time, t9, rho)


def exponential_thermo(t9_0: float = 0.20, rho_0: float = 1.5e4, tau: float = 0.2):
    """Return the exponential history as an analytic ``f(t, y) -> (T9, rho)``.

    The solver is handed the closed-form expressions rather than an interpolated
    table, so the exponential model carries no tabulation error at all.
    """

    def thermo(t: float, y=None):
        t = max(float(t), 0.0)
        return (
            max(t9_0 * np.exp(-t / (3.0 * tau)), T9_FLOOR),
            max(rho_0 * np.exp(-t / tau), RHO_FLOOR),
        )

    return thermo


def write_reference_trajectory(path: Path = REFERENCE_TRAJECTORY, n: int = 4000) -> Path:
    """Generate the reference nova trajectory file described in the docstring."""
    # Logarithmic overall, with a fine linear window across the runaway so that
    # the peak is not softened by the linear interpolation between table rows.
    time = np.unique(
        np.concatenate((
            [0.0],
            np.logspace(-5.0, np.log10(NOVA_T_END), n - 1),
            np.linspace(0.5 * NOVA_T_PEAK, 5.0 * NOVA_T_PEAK, 2000),
        ))
    )

    rising = time <= NOVA_T_PEAK
    fraction = (time[rising] / NOVA_T_PEAK) ** NOVA_RUNAWAY_INDEX
    t9 = np.empty_like(time)
    rho = np.empty_like(time)
    t9[rising] = NOVA_T9_INITIAL + (NOVA_T9_PEAK - NOVA_T9_INITIAL) * fraction
    rho[rising] = NOVA_RHO_INITIAL * (1.0 - 0.5 * fraction)

    rho_peak = 0.5 * NOVA_RHO_INITIAL
    expansion = 1.0 + (time[~rising] - NOVA_T_PEAK) / NOVA_EXPANSION_TIME
    rho[~rising] = rho_peak * expansion ** -3.0
    t9[~rising] = NOVA_T9_PEAK * expansion ** -2.0

    t9 = np.maximum(t9, T9_FLOOR)
    rho = np.maximum(rho, RHO_FLOOR)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("# Reference nova thermodynamic trajectory\n")
        handle.write("# Reconstructed from the description in the proposal; see\n")
        handle.write("# src/thermodynamics.py for the parameterisation.\n")
        handle.write("# time[s]  T9  rho[g/cm^3]\n")
        for t, temperature, density in zip(time, t9, rho):
            handle.write(f"{t:.8e} {temperature:.8e} {density:.8e}\n")
    return path


def nova_trajectory(path: Path = REFERENCE_TRAJECTORY) -> Trajectory:
    """Read the nova trajectory, generating the reference file if it is absent."""
    if not path.exists():
        write_reference_trajectory(path)
    return read_trajectory(path)


def timescales(trajectory: Trajectory) -> tuple[np.ndarray, np.ndarray]:
    """Return the thermodynamic timescales ``|T9 / dT9dt|`` and ``|rho / drhodt|``."""
    time = trajectory.time
    with np.errstate(divide="ignore", invalid="ignore"):
        dt9 = np.gradient(trajectory.t9, time)
        drho = np.gradient(trajectory.rho, time)
        tau_t = np.abs(trajectory.t9 / dt9)
        tau_rho = np.abs(trajectory.rho / drho)
    return np.nan_to_num(tau_t, nan=np.inf), np.nan_to_num(tau_rho, nan=np.inf)
