# Results

Every number below comes from the runs in `results/`, computed with NucNetPy.
Figures are in `figures/`; the LaTeX versions of the tables are in
`results/tables/`.

## 1. The two reference calculations

| Model | Thermodynamic behaviour | $R^{\rm initial}$ | $R^{\rm final}$ | $f_{\rm enh}$ | $t_{\rm fo}$ | $T_{9,\max}$ |
|---|---|---|---|---|---|---|
| Exponential | immediate cooling | 3.68e-3 | **0.460** | 125 | 0.23 s | 0.200 |
| Trajectory | delayed temperature peak | 3.68e-3 | **2.403** | 654 | 152 s | 0.446 |

Both agree with the values quoted in the proposal. The proposal gives
$R^{\rm final} \approx 0.4$–$0.5$ for the exponential model and "of order a
few" for the trajectory model; the calculations give 0.460 and 2.40. The shapes
of the curves also match the proposal's figures: an early bump, a sharp dip, a
long climb, an overshoot, and a flat freeze-out plateau
(`figures/fig07_ratio_comparison.png`).

**One point of bookkeeping.** The methodology chapter defines the ratio from
molar abundances,
$R = (Y_{15\rm N}+Y_{15\rm O})/(Y_{14\rm N}+Y_{14\rm O})$, which for the solar
starting composition gives $R^{\rm initial} = 3.68\times10^{-3}$. Chapter 3
quotes $3.9\times10^{-3}$, which is the same composition evaluated as a ratio of
*mass fractions*, $X(^{15}{\rm N})/X(^{14}{\rm N})$. The two differ by the
factor $14/15$. The molar definition is used throughout here, because that is
the one the methodology chapter gives; the enhancement factor is barely
affected either way.

## 2. Why the ratio moves when it does

The reaction flows (`fig08`, `fig09`) and the timescales (`fig10`, `fig11`)
explain the shape of the curves.

**The early dip is real, not numerical.** At the starting temperature the two
fastest reactions in the whole network are $^{18}{\rm O}({\rm p},\alpha)^{15}{\rm N}$
and $^{15}{\rm N}({\rm p},\alpha)^{12}{\rm C}$. The second destroys the solar
$^{15}$N far faster than anything replaces it, so $R$ falls by more than an
order of magnitude before the CNO cycle has begun to turn. This is also what
sets the initial nuclear energy generation, about $10^{15}$ erg g$^{-1}$
s$^{-1}$ (`fig14`, `fig15`).

**The climb is the hot-CNO cycle turning over.** $^{12}{\rm C}({\rm p},\gamma)$
feeds $^{13}$N, and from there the flow reaches $^{14}$O and $^{15}$O. Because
$^{14}$O and $^{15}$O are beta-unstable with half-lives of 71 s and 122 s, they
build up while the material is hot and then decay to $^{14}$N and $^{15}$N once
it is not. $^{15}$O is made faster than $^{14}$O, so the ratio ends far above
its starting value.

**Freeze-out is when the proton-capture timescale overtakes the
thermodynamic timescale.** In the exponential model
$\tau_T = 3\tau = 0.6$ s throughout, and the proton-capture timescales rise
through it within the first second; the ratio is within one per cent of its
final value by $t_{\rm fo} = 0.23$ s. In the trajectory model $\tau_T$ is
effectively infinite during the quiescent phase, collapses to a few seconds
across the runaway, and then grows again; freeze-out follows at
$t_{\rm fo} = 152$ s, once the beta decays of $^{14}$O and $^{15}$O have run to
completion.

**The trajectory produces more $^{15}$N because it burns hotter, for longer.**
The peak temperature is 0.446 against 0.200, and material sits near the peak
for tens of seconds instead of cooling away from it immediately.

## 3. Steady flow

The flow ratios $Q_{ab} = F_a / F_{^{14}{\rm N}({\rm p},\gamma)}$ (`fig12`,
`fig13`) answer the quasi-equilibrium question directly: **the CNO cycle does
not reach steady flow, except momentarily.**

In the trajectory model the three ratios span more than six orders of magnitude
during the quiescent phase, converge towards unity only in the few seconds
around the temperature peak, and separate again immediately afterwards. So
there is a brief interval near maximum temperature where the cycle approaches
quasi-steady flow, and nothing resembling it before or after.
$^{14}{\rm N}({\rm p},\gamma)^{15}{\rm O}$ is the slowest link throughout,
which is the classical CNO bottleneck.

## 4. Sensitivity to the expansion timescale

| $\tau$ (s) | $R^{\rm final}$ | $f_{\rm enh}$ | $t_{\rm fo}$ (s) |
|---|---|---|---|
| 0.05 | 0.584 | 159 | 0.052 |
| 0.10 | 0.609 | 166 | 0.112 |
| 0.20 | 0.460 | 125 | 0.226 |
| 0.50 | 0.330 | 90 | 0.454 |
| 1.00 | 0.315 | 86 | 0.352 |

The final ratio is not monotonic in $\tau$ (`fig17`). It peaks near
$\tau \approx 0.1$ s and falls away on both sides. Fast cooling freezes the
composition before much $^{15}$O can be made; slow cooling leaves enough time
at moderate temperature for $^{15}{\rm N}({\rm p},\alpha)$ to destroy some of
what was made. The largest enhancement comes from the intermediate case, and
the whole range is only a factor of two wide — much smaller than the factor of
five between the exponential and trajectory models.

## 5. Network size

See `results/summary.json` and `results/tables/tab_network_size.tex` for the
recorded values, and `figures/fig16_network_size.png` for the curves.

## 6. Numerical quality

| Check | Result |
|---|---|
| Baryon number, every reaction | exact, all three networks |
| Charge changes | only $\pm1$, i.e. beta decays with the leptons left out of the ReacLib record |
| Dead-end nuclides | none |
| Infinite or prompt rates in range | none |
| Beta-decay rates vs measured half-lives | within 0.2% for $^{13}$N, $^{14}$O, $^{15}$O, $^{17}$F, $^{18}$F |
| $\sum_i A_i Y_i$ at the end of a run | 1 to within 1e-15 (exponential), 1e-15 (trajectory) |
| Tolerance convergence | $R^{\rm final}$ stable to six figures between `rtol` 1e-6 and 1e-8 |

Three repairs were needed to the network before it behaved:

1. **Particle-unbound nuclides were dead ends.** ReacLib produces $^{5}$Li,
   $^{8}$Be and $^{9}$B but gives them no way to decay, so they slowly
   accumulated material that should have gone back to $\alpha$ particles and
   protons. Every reaction producing one is now written directly to its
   break-up products, the way $^{8}$B beta decay is normally written as
   $^{8}{\rm B} \rightarrow 2\alpha$.

2. **One prompt proton emitter was left in.** $^{18}{\rm Na} \rightarrow
   {\rm p} + ^{17}{\rm Ne}$ carries a rate of $5\times10^{20}$ s$^{-1}$ in
   ReacLib. It was eliminated the same way. This is found automatically now,
   by testing for one-body rates above $10^{12}$ s$^{-1}$.

3. **The ReacLib fits must not be extrapolated below $T_9 = 0.01$.** Below it
   some neutron-induced rates run away to $10^{24}$ s$^{-1}$ and beyond. Both
   thermodynamic histories are floored at $T_9 = 0.01$, which costs nothing:
   charged-particle reactions are already dead there, and the beta decays that
   still matter during the cooling tail do not depend on temperature.

None of the three changed the diagnostic ratio — $R^{\rm final}$ for the
exponential model was 0.46004 before and after — but each was a defect that
would have been fair to ask about.

## 7. What is a modelling assumption, not a result

- **The nova trajectory is a reconstruction.** No hydrodynamic trajectory file
  was available, so `data/trajectories/nova_reference.txt` was built to match
  the description in the proposal: quiescence at $T_9 = 0.091$ and
  $\rho = 2.21\times10^4$ g cm$^{-3}$, a runaway peaking at $T_9 = 0.447$ near
  100 s, then adiabatic expansion. The trajectory-model numbers are therefore
  conditional on that history. The file is plain three-column text and the
  analysis reads it through `nucnetpy.read_trajectory`, so substituting a real
  one requires no code changes.
- **Deuterium is set to zero** in the initial composition. Material accreted
  from the companion has no primordial deuterium left; including the solar
  value puts a spurious spike on the energy generation in the first
  microsecond and changes nothing else.
- **No screening.** At $T_9 \leq 0.45$ and $\rho \leq 2.2\times10^4$ g
  cm$^{-3}$ the plasma is weakly coupled and screening changes the CNO rates by
  well under one per cent. NucNetPy's `SkyNetScreening` can be switched on if a
  quantitative bound is wanted.
- **Solar starting composition**, as the proposal specifies. Real nova envelopes
  are mixed with white-dwarf material, which would raise the CNO abundances by
  a large factor and change the energetics, though not obviously the ratio.
