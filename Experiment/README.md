# Nova hydrogen-burning experiment, run with NucNetPy

This directory contains the calculations described in the proposal, carried out
with [NucNetPy](https://github.com/mtalafha90/NucNetPy), the pure-Python nuclear
reaction-network package. Nothing in the NucNetPy repository was changed: all
the code here reads from its public interface.

## What is calculated

The computational plan of the methodology chapter, in order:

1. the reference exponential single-zone calculation,
   `T9_0 = 0.20`, `rho_0 = 1.5e4 g/cm^3`, `tau = 0.2 s`, followed to 100 s;
2. the reference trajectory-based single-zone calculation, followed to
   `3.15e7 s`;
3. the comparison of the two histories through the diagnostic ratio
   `R = (Y(15N) + Y(15O)) / (Y(14N) + Y(14O))`;
4. the trajectory calculation repeated on three networks (`Z<=10`, `Z<=20`,
   `Z<=30`) to test whether the result depends on network size;
5. the exponential calculation repeated for
   `tau = 0.05, 0.10, 0.20, 0.50, 1.00 s`;
6. the hot-CNO reaction flows, the nuclear and thermodynamic timescales, the
   freeze-out times, and the flow ratios that test for quasi-steady flow.

## Layout

```
src/
  build_networks.py   ReacLib + mass evaluation -> the network archives (run once)
  network_io.py       archive -> nucnetpy.Network
  composition.py      solar initial composition -> molar abundances
  thermodynamics.py   the exponential and trajectory temperature-density histories
  diagnostics.py      R_15/14, freeze-out, flows, timescales, steady-flow ratios
  run_experiment.py   the driver: every run, every output file
  make_figures.py     figures and LaTeX tables from the stored results
data/
  networks/           the three networks, as JSON archives and NucNetPy XML
  trajectories/       the nova temperature-density history
results/
  *_evolution.csv     T9, rho, R_15/14 and the tracked isotopes for every run
  *_flows.csv         reaction flows, timescales and energy generation
  summary.json        every diagnostic number, per run
  tables/*.tex        LaTeX tables ready to be included in the proposal
figures/              the figures
RESULTS.md            what the calculations show
```

## Reproducing

```bash
python -m pip install nucnetpy matplotlib      # numpy and scipy come with nucnetpy
python src/run_experiment.py                   # all runs; a few hours on one core
python src/make_figures.py                     # figures and tables
```

Individual runs can be given by name, and are independent, so they can be run in
parallel:

```bash
python src/run_experiment.py exp_ref
python src/run_experiment.py traj_z30
```

Rebuilding the networks additionally needs `pynucastro`, which supplies the JINA
ReacLib snapshot and the atomic mass evaluation:

```bash
python -m pip install pynucastro
python src/build_networks.py
```

This step is not needed to repeat the calculations: the network archives in
`data/networks/` are complete and are read directly.

## Physics inputs

| Ingredient | Choice |
|---|---|
| Rates | JINA ReacLib snapshot `20180319default2`, seven-parameter fits, forward and reverse |
| Masses | AME atomic mass evaluation |
| Composition | Solar, Bergemann, Lodders & Palme (2025), with solar isotopic splits |
| Nuclides | `Z <= Z_max` and neutron excess `N - Z <= 3`; drip-line nuclides without a measured mass are dropped |
| Solver | `nucnetpy.evolve_zone`, SciPy BDF with NucNetPy's analytic Jacobian, `rtol = 1e-8`, `atol = 1e-30` |
| Screening | None. At `T9 <= 0.45` and `rho <= 2e4 g/cm^3` the plasma is weakly coupled and screening changes the CNO rates by well under a per cent |

## A note on the trajectory file

The nova temperature-density history is `data/trajectories/nova_reference.txt`.
It is a **reconstruction**, built to match the behaviour described in the
proposal (a quiescent phase at `T9 = 0.091` and `rho = 2.21e4 g/cm^3`, a
runaway peaking at `T9 = 0.447` about 100 s in, then adiabatic expansion), not
an extract from a hydrodynamic nova model. It is kept as a separate data file
in a plain three-column `time T9 rho` format, which is what
`nucnetpy.read_trajectory` reads, so a trajectory from a real nova model can be
dropped in its place and everything downstream will run unchanged. The
parameterisation is documented in `src/thermodynamics.py`.
