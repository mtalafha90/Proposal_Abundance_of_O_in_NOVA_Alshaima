"""
Generate a literature-constrained synthetic thermodynamic trajectory
inspired by Iliadis et al. (2002) model S1.

Only the model identity and published peak temperature are taken from
the paper. The time and density evolution are explicit analytic
assumptions and are not the original hydrodynamic trajectory.
"""

from pathlib import Path
import numpy as np

MODEL_TPEAK_GK = 0.418

T_PEAK_S = 100.0
T0_GK = 0.070
T_FLOOR_GK = 0.010
RHO0 = 2.20e4
RHO_AT_TPEAK = 4.00e3
TAU_T = 25.0
ALPHA_T = 1.25
TAU_RHO = 45.0
ALPHA_RHO = 0.90
LOGISTIC_MIDPOINT = 0.86
LOGISTIC_STEEPNESS = 16.0
FINAL_TIME_S = 3000.0


def normalised_logistic(u, midpoint, steepness):
    raw = 1.0 / (1.0 + np.exp(-steepness * (u - midpoint)))
    raw0 = 1.0 / (1.0 + np.exp(-steepness * (0.0 - midpoint)))
    raw1 = 1.0 / (1.0 + np.exp(-steepness * (1.0 - midpoint)))
    return (raw - raw0) / (raw1 - raw0)


def build_time_grid():
    a = np.arange(0.0, 60.0, 0.5)
    b = np.arange(60.0, 140.0, 0.05)
    c = np.arange(140.0, 300.0, 0.25)
    d = np.geomspace(300.0, FINAL_TIME_S, 260)
    return np.unique(
        np.concatenate([a, b, c, d, np.array([T_PEAK_S, FINAL_TIME_S])])
    )


def build_profile(t):
    T9 = np.empty_like(t)
    rho = np.empty_like(t)

    heating = t <= T_PEAK_S
    cooling = ~heating

    u = np.clip(t[heating] / T_PEAK_S, 0.0, 1.0)
    s = normalised_logistic(u, LOGISTIC_MIDPOINT, LOGISTIC_STEEPNESS)

    T9[heating] = T0_GK + (MODEL_TPEAK_GK - T0_GK) * s
    rho[heating] = np.exp(
        np.log(RHO0) + (np.log(RHO_AT_TPEAK) - np.log(RHO0)) * s
    )

    dt = t[cooling] - T_PEAK_S
    T9[cooling] = T_FLOOR_GK + (MODEL_TPEAK_GK - T_FLOOR_GK) * np.exp(
        -((dt / TAU_T) ** ALPHA_T)
    )
    rho[cooling] = RHO_AT_TPEAK * np.exp(
        -((dt / TAU_RHO) ** ALPHA_RHO)
    )
    rho[cooling] = np.maximum(rho[cooling], 1.0e-12)

    ipeak = np.argmin(np.abs(t - T_PEAK_S))
    T9[ipeak] = MODEL_TPEAK_GK
    rho[ipeak] = RHO_AT_TPEAK

    return T9, rho


def main():
    t = build_time_grid()
    T9, rho = build_profile(t)

    out = Path("iliadis2002_S1_synthetic_benchmark.txt")
    header = """Literature-constrained synthetic S1-like nova trajectory
NOT the original hydrodynamic S1 time series.
Columns: time_s temperature_T9_GK density_g_cm3"""
    np.savetxt(
        out,
        np.column_stack([t, T9, rho]),
        fmt="%.8e",
        header=header,
    )
    print(f"Wrote {out} with {len(t)} rows")


if __name__ == "__main__":
    main()
