# Nova hydrogen-burning experiment, run with NucNetPy

This directory contains the calculations described in the proposal, carried out
with [NucNetPy](https://github.com/mtalafha90/NucNetPy), the pure-Python nuclear
reaction-network package. Nothing in the NucNetPy repository was changed: all
the code here reads from its public interface.

NucNetPy is a reimplementation, in Python, of the approach taken by the C/C++
[NucNet Tools](https://sourceforge.net/projects/nucnet-tools/) of B. S. Meyer
and the Webnucleo group at Clemson University, from which it inherits its data
model and its JINA/libnucnet-compatible input formats. The earlier solar
abundance studies from this group used NucNet Tools itself; this study is the
first application of NucNetPy.

## What is calculated

The computational plan of the methodology chapter, in order:

1. the reference exponential single-zone calculation,
   `T9_0 = 0.20`, `rho_0 = 1.5e4 g/cm^3`, `tau = 0.2 s`, followed to 100 s;
2. the reference trajectory-based single-zone calculation, following the
   nova profile to the end of its file at `3000 s`;
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
  validate_networks.py checks on the networks and the composition, before any run
  run_experiment.py   the driver: every run, every output file
  make_figures.py     figures and LaTeX tables from the stored results
data/
  networks/           the three networks, as JSON archives and NucNetPy XML
  trajectories/       the nova temperature-density history
results/
  *_evolution.csv     T9, rho, R_15/14 and the tracked isotopes for every run
  *_flows.csv         reaction flows, timescales and energy generation
  summary.json        every diagnostic number, per run
  network_validation.json  the checks, recorded
  tables/*.tex        LaTeX tables ready to be included in the proposal
figures/              the figures
RESULTS.md            what the calculations show
```

## Reproducing

```bash
python -m pip install nucnetpy matplotlib      # numpy and scipy come with nucnetpy
python src/validate_networks.py                # check the networks first
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
| Nuclides | `Z <= Z_max` and neutron excess `N - Z <= 3 + Z/4`, a band wide enough to keep every stable isotope of the starting composition; drip-line nuclides without a measured mass, and particle-unbound nuclides, are eliminated |
| Solver | `nucnetpy.evolve_zone`, SciPy BDF with NucNetPy's analytic Jacobian, `rtol = 1e-8`, `atol = 1e-25` |
| Screening | None. At `T9 <= 0.45` and `rho <= 2e4 g/cm^3` the plasma is weakly coupled and screening changes the CNO rates by well under a per cent |

## The trajectory file

The nova temperature-density history is
`data/trajectories/iliadis2002_S1_synthetic_benchmark.txt`: three columns of
time (s), temperature (in units of `1e9 K`) and density (g cm^-3), read by
`nucnetpy.read_trajectory`. It is a literature-constrained synthetic benchmark
anchored to the peak temperature of model S1 of Iliadis et al. (2002), not a
hydrodynamic trajectory; the generator and its parameters are in
`data/trajectories/provenance/`.

```
t = 0          T9 = 0.070     rho = 2.200e4 g/cm^3
t = 100.0 s    T9 = 0.418     rho = 4.000e3        <- hottest point
t = 3000 s     T9 = 0.010     rho = 1.0e-12        <- last row
```

The burning episode is a broad peak rather than a spike: the temperature stays
above `T9 = 0.2` for 38.4 s and above `T9 = 0.1` for 64.2 s.

Two things are worth knowing about how it is used.

- **The runs stop where the file stops**, at `t = 3000 s`, rather than at the
  `3.15e7 s` named in the proposal. Going further would mean holding the last
  row's temperature and density fixed for four more decades of time, which is
  an extrapolation rather than a result. Nothing is lost: the diagnostic ratio
  is within one per cent of its final value by 144.5 s.
- **The temperature is floored at `T9 = 0.01`** because ReacLib's fits are not
  made for temperatures below that and some of them diverge if pushed there.
  For this profile the floor never binds: the history is generated with the
  same floor and approaches it asymptotically from above, so no tabulated value
  is altered.
