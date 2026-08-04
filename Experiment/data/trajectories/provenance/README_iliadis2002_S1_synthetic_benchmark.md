# Iliadis et al. (2002) S1-like synthetic nova trajectory

## Status

This is a **literature-constrained synthetic benchmark**, not the
unpublished original hydrodynamic S1 temperature-density-time profile.

The paper supplies the S1 model identity and global parameters but does
not tabulate the time series. The authentic source-derived anchors used
here are:

- model: S1
- white dwarf: ONe
- white-dwarf mass: 1.35 solar masses
- assumed core-envelope mixing: 50 per cent
- peak temperature: 0.418 GK
- accretion rate: 1.6e-10 solar masses per year

## Analytic assumptions

The following quantities are modelling choices made only to produce a
usable NucNetPy benchmark:

- peak time: 100.0 s
- initial temperature: 0.070 GK
- initial density: 2.200e+04 g cm^-3
- density at maximum temperature:
  4.000e+03 g cm^-3
- stretched-exponential cooling and expansion after the maximum
- final tabulated time: 3000 s

These assumptions must not be cited as results of Iliadis et al. (2002).

## File format

`iliadis2002_S1_synthetic_benchmark.txt` contains:

1. time in seconds
2. temperature T9 in GK
3. density in g cm^-3

The file has 2621 rows and SHA-256:

`420ad5726205b9e949b7e33a4cef32f0da6f4dea1f31d6f3cb0cde41fa793c72`

## Derived profile properties

- peak temperature: 0.418000 GK
- time of maximum temperature: 100.000 s
- density at maximum temperature:
  4.000000e+03 g cm^-3
- time above T9=0.1: 64.15 s
- time above T9=0.2: 38.35 s
- time above T9=0.3: 22.05 s

## Use in a paper

Describe it as an analytic or synthetic S1-like benchmark anchored to
the published peak temperature. Do not call it the Iliadis S1
hydrodynamic trajectory.

## References

C. Iliadis et al., 2002, The Astrophysical Journal Supplement Series,
142, 105-137, DOI 10.1086/341400.
