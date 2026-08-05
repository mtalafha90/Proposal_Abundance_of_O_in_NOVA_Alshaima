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
The nova trajectory is
``data/trajectories/iliadis2002_S1_synthetic_benchmark.txt``, read with
:func:`nucnetpy.read_trajectory`.  It starts at ``T9 = 0.070`` and
``rho = 2.200e4 g/cm^3``, rises to ``T9 = 0.418`` at ``t = 100.0 s``, where the
density is ``4.000e3 g/cm^3``, and then cools and expands out to
``t = 3000 s``.  The temperature stays above ``T9 = 0.2`` for 38.35 s and above
``T9 = 0.1`` for 64.15 s, so the burning episode is a broad peak rather than a
spike.

What this file is, and is not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is a *literature-constrained synthetic benchmark*, not a hydrodynamic
trajectory.  Only the global parameters are taken from the published S1 model
of Iliadis et al. (2002, ApJS 142, 105): an ONe white dwarf of 1.35 solar
masses, 50 per cent core-envelope mixing, and a peak temperature of 0.418 GK.
Iliadis et al. do not tabulate the time series, so the *shape* of this profile
-- the 100 s rise, the initial and peak densities, and the stretched-exponential
cooling and expansion that follow -- is an analytic construction.  Those
choices are documented in ``data/trajectories/provenance/``, which holds the
generator script, its parameter summary and the upstream README.

The consequence for this study is that both thermodynamic histories are now
analytic.  The trajectory model is a more elaborate parameterisation than the
exponential one -- it has a heating phase, a broad maximum and different
cooling laws for temperature and density -- but it is not evidence of what a
hydrodynamic calculation would do, and no result here should be described as
one.  The published anchor is the peak temperature alone.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from nucnetpy import Trajectory, read_trajectory

ROOT = Path(__file__).resolve().parent.parent
TRAJECTORY_DIR = ROOT / "data" / "trajectories"

#: The nova profile in use: time (s), T9, density (g/cm^3).  See the module
#: docstring -- this is a synthetic benchmark anchored to the published peak
#: temperature of the Iliadis et al. (2002) S1 model, not a hydrodynamic run.
NOVA_PROFILE = TRAJECTORY_DIR / "iliadis2002_S1_synthetic_benchmark.txt"


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



def nova_trajectory(path: Path | None = None) -> Trajectory:
    """Read the nova trajectory, floored so the rate fits stay in range.

    :data:`NOVA_PROFILE` is read unless an explicit ``path`` is given.

    The floors applied here are a safeguard, and for the S1 benchmark they never
    bind.  That profile is generated with the same ``T9 = 0.01`` floor and
    approaches it asymptotically from above -- within 1e-5 of it by
    ``t = 250 s``, and equal to it in the stored file, which carries eight
    significant figures, from ``t = 406 s``.  Its lowest density, ``1e-12
    g/cm^3``, stays eighteen orders of magnitude above :data:`RHO_FLOOR`.
    Clamping therefore changes nothing here, as
    ``max(profile - floored profile)`` confirms: it is zero for both columns.

    Do not confuse :data:`RHO_FLOOR` with the ``1e-12 g/cm^3`` lower bound that
    the generator applies while *building* the profile.  That one does bind, for
    the last 21 rows from ``t = 2511 s`` onwards, but it is part of the
    definition of the trajectory rather than something imposed here, and it
    takes effect long after the composition has frozen at ``t = 144.5 s``.

    The floors would matter for a profile that really did run below the range
    the ReacLib fits were made for, where some of them diverge.  In that case
    nothing physical is lost, because charged-particle reactions are long dead
    at ``T9 = 0.01`` and the beta decays that still run do not depend on
    temperature at all.
    """
    source = path if path is not None else NOVA_PROFILE
    if not source.exists():
        raise FileNotFoundError(
            f"nova profile not found: {source}.  The trajectory is data, not "
            "something this module can synthesise; supply the file or pass an "
            "explicit path."
        )
    trajectory = read_trajectory(source)

    # Drop repeated time rows.  The S1 benchmark file carries the temperature
    # maximum twice, at t = 100 s, because its generator appends the peak time
    # to a grid that already contains it and the two differ only below the
    # precision the file is written with.  Interpolation is untroubled by that,
    # and nothing in the production pipeline differentiates the table -- the
    # thermodynamic timescales are taken on the output grid -- but a zero gap
    # would give an infinite derivative to anything that did, including
    # :func:`timescales` below.  Removing the duplicate leaves the interpolant
    # unchanged.
    keep = np.concatenate(([True], np.diff(trajectory.time) > 0.0))

    return Trajectory(
        trajectory.time[keep],
        np.maximum(trajectory.t9[keep], T9_FLOOR),
        np.maximum(trajectory.rho[keep], RHO_FLOOR),
    )


def peak_time(trajectory: Trajectory) -> float:
    """Time at which the trajectory is hottest."""
    return float(trajectory.time[int(np.argmax(trajectory.t9))])


def timescales(trajectory: Trajectory) -> tuple[np.ndarray, np.ndarray]:
    """Return the thermodynamic timescales ``|T9 / dT9dt|`` and ``|rho / drhodt|``."""
    time = trajectory.time
    with np.errstate(divide="ignore", invalid="ignore"):
        dt9 = np.gradient(trajectory.t9, time)
        drho = np.gradient(trajectory.rho, time)
        tau_t = np.abs(trajectory.t9 / dt9)
        tau_rho = np.abs(trajectory.rho / drho)
    return np.nan_to_num(tau_t, nan=np.inf), np.nan_to_num(tau_rho, nan=np.inf)
